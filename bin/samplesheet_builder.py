#!/bin/env python3
import os
import csv
from collections import defaultdict
import sys
import re
import argparse

# Set default directories
DEFAULT_SAMPLESHEET_DIR = "/data/khanlab/projects/DATA/Sequencing_Tracking_Master"
DEFAULT_INPUT_DIR = "/data/khanlab/projects/DATA"
DEFAULT_BAM_DIR = "/data/khanlab3/David_Milewski_StJude_202504"

# Check required positional args (allow optional flags like --genome)
if len(sys.argv) < 3:
    print(
        f"Usage: python {sys.argv[0]} <patient_id> <case_name> [--genome <hg19|hg38|mm39>]"
    )
    print(f"Default Mastersheet Directory: {DEFAULT_SAMPLESHEET_DIR}")
    print(f"Default Input Directory: {DEFAULT_INPUT_DIR}")
    print("To use custom directories, modify the script:")
    print(f"   - Change 'DEFAULT_SAMPLESHEET_DIR' to your samplesheet directory path")
    print(f"   - Change 'DEFAULT_INPUT_DIR' to your input directory path")
    sys.exit(1)

# Extract required positionals
sample_id = sys.argv[1]
case_name = sys.argv[2]

# Optional: --genome fallback (minimal parser, no argparse wiring)
fallback_genome = None
if "--genome" in sys.argv:
    i = sys.argv.index("--genome")
    if i + 1 < len(sys.argv):
        fallback_genome = (sys.argv[i + 1] or "").strip().lower()
        # normalize a few common aliases
        if fallback_genome in {"grch38", "38"}:
            fallback_genome = "hg38"
        elif fallback_genome in {"grch37", "b37", "hs37d5", "37", "19"}:
            fallback_genome = "hg19"
        elif fallback_genome in {"grcm39"}:
            fallback_genome = "mm39"

# Use default directories
samplesheet_dir = DEFAULT_SAMPLESHEET_DIR
inputdir = DEFAULT_INPUT_DIR

print(f"Sample ID: {sample_id}")
print(f"Case Name: {case_name}")
print(f"Samplesheet Directory: {samplesheet_dir}")
print(f"Input Directory: {inputdir}")
print(f"BAM Directory: {DEFAULT_BAM_DIR}")
if fallback_genome:
    print(f"Fallback genome from CLI: {fallback_genome}")


def read_and_map_samplesheet(
    samplesheet, inputdir, column_mapping, sample_id, case_name
):
    samplesheet_data = []
    invalid_paths = []
    try:
        with open(
            samplesheet, "r", newline="", encoding="utf-8", errors="replace"
        ) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                # column mapping
                mapped_row = {
                    new_column: row.get(old_column, "").strip()
                    for old_column, new_column in column_mapping.items()
                }
                if "type" in mapped_row:
                    mapped_row["type"] = mapped_row["type"].replace(" ", "_")
                if "type" in mapped_row and mapped_row["type"] == "blood_DNA":
                    mapped_row["type"] = "normal_DNA"
                if "Diagnosis" in mapped_row:
                    mapped_row["Diagnosis"] = re.sub(
                        r",\s*", "_", mapped_row["Diagnosis"]
                    )
                    mapped_row["Diagnosis"] = mapped_row["Diagnosis"].replace(" ", ".")

                # --- expand multiple casenames (append once) ---
                appended = False
                if "," in mapped_row.get("casename", ""):
                    casenames = mapped_row["casename"].split(",")
                    for individual_casename in casenames:
                        new_row = dict(mapped_row)
                        new_row["casename"] = individual_casename.strip()
                        samplesheet_data.append(new_row)
                    appended = True
                else:
                    pass

                # --- expand multi seq_type (append once) ---
                if "/" in mapped_row.get("seq_type", ""):
                    # add the primary as-is on first append
                    if not appended:
                        samplesheet_data.append(mapped_row)
                        appended = True
                    seqtypes = mapped_row["seq_type"].split("/")
                    # start from the second item (index 1)
                    for individual_seqtype in seqtypes[1:]:
                        new_row = dict(mapped_row)
                        new_row["seq_type"] = individual_seqtype.strip()
                        samplesheet_data.append(new_row)

                # if nothing special, append once
                if not appended:
                    samplesheet_data.append(mapped_row)

    except UnicodeDecodeError:
        print(f"Error decoding file: {samplesheet}. Please check the file encoding.")
        return [], []

    ALLOWED_SAMPLE_CAPTURES = {
        "access",
        "polya_stranded",
        "polya",
        "ribozero",
        "smartrna",
        "ribodepleted_nebnext_v2",
        "clin.ex.v1",
        "seqcapez.hu.ex.v3",
        "seqcapez.rms.v1",
        "agilent.v7",
        "idt_v2_plus",
        "xgen-hyb-panelv2",
        "comp_ex_v1",
        "seqcapez.hu.ex.utr.v1",
    }

    # Filter rows matching sample_id and case_name
    filtered_samplesheet_data = []
    for row in samplesheet_data:
        if (
            row.get("sample") == sample_id
            and row.get("casename") == case_name
            and row.get("seq_type") in ["E-il", "P-il", "T-il"]
        ):
            if row.get("type") == "normal_RNA":
                row["type"] = "cell_line_RNA"

            sample_captures = row.get("sample_captures", "").strip().lower()
            row["sample_captures"] = sample_captures
            if not sample_captures:
                print(
                    f"ERROR: 'sample_captures' is missing or empty for sample {sample_id}, case {case_name}"
                )
                sys.exit(1)
            if sample_captures not in ALLOWED_SAMPLE_CAPTURES:
                print(
                    f"ERROR: Unrecognized 'sample_captures' value '{sample_captures}' for sample {sample_id}, case {case_name}."
                )
                print(
                    "Please add this sample capture type to the workflow before creating the samplesheet."
                )
                sys.exit(1)

            # default genome if missing or invalid; prefer CLI fallback when provided

            if fallback_genome:
                row["genome"] = fallback_genome
                print(f"genome missing in sheet, using CLI fallback: {fallback_genome}")
            else:
                if row.get("genome") not in ["hg19", "hg38", "mm39"]:
                    print("genome information missing, defaulting to hg19")
                    row["genome"] = "hg19"

            library_id = row["library"]
            if "." in library_id:
                print(
                    f"ERROR: Library ID '{library_id}' contains a '.' character. "
                    "This will cause issues with visualizing results on the website.\n"
                    "Please rename the library ID to replace '.' with '_' and relaunch the script."
                )
                sys.exit(1)
            fcid = row.get("FCID", "")

            # Compose FASTQ expectations
            if fcid:
                read1 = f"{inputdir}/Sample_{library_id}_{fcid}/Sample_{library_id}_{fcid}_R1.fastq.gz"
                read2 = f"{inputdir}/Sample_{library_id}_{fcid}/Sample_{library_id}_{fcid}_R2.fastq.gz"
            else:
                read1 = (
                    f"{inputdir}/Sample_{library_id}/Sample_{library_id}_R1.fastq.gz"
                )
                read2 = (
                    f"{inputdir}/Sample_{library_id}/Sample_{library_id}_R2.fastq.gz"
                )

            # Check FASTQs; else fall back to BAM
            if os.path.exists(read1) and os.path.exists(read2):
                row["read1"] = read1
                row["read2"] = read2
                row["bam"] = ""  # ensure bam column present when FASTQs are used
                filtered_samplesheet_data.append(row)
            else:
                bam_path = os.path.join(DEFAULT_BAM_DIR, f"{library_id}.bam")
                if os.path.exists(bam_path):
                    row["read1"] = ""
                    row["read2"] = ""
                    row["bam"] = bam_path
                    filtered_samplesheet_data.append(row)
                else:
                    # Try alternate BAM naming conventions
                    alt_library_id = library_id.replace("_Exome", ".Exome").replace(
                        "_RNA-Seq", ".RNA-Seq"
                    )
                    alt_bam_path = os.path.join(
                        DEFAULT_BAM_DIR, f"{alt_library_id}.bam"
                    )

                    if os.path.exists(alt_bam_path):
                        row["read1"] = ""
                        row["read2"] = ""
                        row["bam"] = alt_bam_path
                        filtered_samplesheet_data.append(row)
                    else:
                        both_bams = f"{bam_path} || {alt_bam_path}"
                        invalid_paths.append((read1, read2, both_bams))

    return filtered_samplesheet_data, invalid_paths


def fill_matches_for_group(rows):
    """
    Strict matching:
      - If BOTH Matched_RNA and Matched_normal are already populated for a DNA row -> skip it.
      - Else fill only missing fields (don't overwrite non-empty):
          * Matched_normal <- normal_DNA library (if exactly one exists)
          * Matched_RNA    <- xeno_RNA or tumor_RNA library (if exactly one total exists)
      - If multiple normal_DNA libs exist (same sample+casename) and any DNA row needs a normal -> exit with error.
      - If multiple RNA libs exist (same sample+casename) and any DNA row needs RNA -> exit with error.
    Assumes library ID is the unique identifier per row.
    """

    def norm(v):
        t = (v or "").strip()
        t = t.replace("-", "_").replace(" ", "_")
        return t.lower()

    # Bucket by type
    by_type = {
        "normal_dna": [],
        "xeno_rna": [],
        "tumor_rna": [],
        "xeno_dna": [],
        "tumor_dna": [],
    }
    for r in rows:
        t = norm(r.get("type"))
        if t in by_type:
            by_type[t].append(r)

    # Unique libs per type
    def libs(rs):
        return sorted(
            {
                (r.get("library") or "").strip()
                for r in rs
                if (r.get("library") or "").strip()
            }
        )

    normal_dna_libs = libs(by_type["normal_dna"])
    xeno_rna_libs = libs(by_type["xeno_rna"])
    tumor_rna_libs = libs(by_type["tumor_rna"])
    all_rna_libs = sorted(set(xeno_rna_libs) | set(tumor_rna_libs))

    def dna_row_needs(rna=False, normal=False):
        for r in rows:
            t = norm(r.get("type"))
            if t not in {"xeno_dna", "tumor_dna"}:
                continue
            mrna = (r.get("Matched_RNA") or "").strip()
            mnorm = (r.get("Matched_normal") or "").strip()
            # Skip rows already fully matched
            if mrna and mnorm:
                continue
            if rna and not mrna:
                return True
            if normal and not mnorm:
                return True
        return False

    need_rna = dna_row_needs(rna=True)
    need_norm = dna_row_needs(normal=True)

    # Strict checks
    errs = []
    sample = (rows[0].get("sample") or "").strip()
    case = (rows[0].get("casename") or "").strip()

    if need_norm and len(normal_dna_libs) > 1:
        errs.append(f"Multiple normal_DNA libraries: {', '.join(normal_dna_libs)}")

    if need_rna and len(all_rna_libs) > 1:
        errs.append(f"Multiple RNA libraries: {', '.join(all_rna_libs)}")

    if errs:
        print(
            f"ERROR: Ambiguous matches for sample '{sample}', case '{case}'. "
            + " | ".join(errs)
            + ". Unable to set matches deterministically. Please select a single library."
        )
        sys.exit(1)

    normal_dna_lib = normal_dna_libs[0] if len(normal_dna_libs) == 1 else ""
    rna_lib = all_rna_libs[0] if len(all_rna_libs) == 1 else ""

    # Fill per DNA row
    for r in rows:
        r_type = norm(r.get("type"))
        if r_type not in {"xeno_dna", "tumor_dna"}:
            continue

        matched_rna = (r.get("Matched_RNA") or "").strip()
        matched_normal = (r.get("Matched_normal") or "").strip()

        # Treat "N/A", "NA", empty strings as missing values
        if matched_rna.upper() in {"N/A", "NA", ""}:
            matched_rna = ""
        if matched_normal.upper() in {"N/A", "NA", ""}:
            matched_normal = ""

        if matched_rna and matched_normal:
            continue

        if not matched_normal and normal_dna_lib:
            r["Matched_normal"] = normal_dna_lib

        if not matched_rna and rna_lib:
            r["Matched_RNA"] = rna_lib


def process_samplesheets(
    samplesheet_dir, inputdir, column_mapping, sample_id, case_name
):
    all_filtered_data = []
    all_invalid_paths = []

    for filename in os.listdir(samplesheet_dir):
        filepath = os.path.join(samplesheet_dir, filename)
        if os.path.isfile(filepath):
            print(f"Processing file: {filepath}")
            filtered_data, invalid_paths = read_and_map_samplesheet(
                filepath, inputdir, column_mapping, sample_id, case_name
            )
            all_filtered_data.extend(filtered_data)
            all_invalid_paths.extend(invalid_paths)

    # Group data and write to output files
    grouped_data = defaultdict(list)
    for row in all_filtered_data:
        grouped_data[(row["sample"], row["casename"])].append(row)

    for key in grouped_data:
        # de-dup per group
        grouped_data[key] = list(
            {tuple(sorted(d.items())): d for d in grouped_data[key]}.values()
        )
        fill_matches_for_group(grouped_data[key])

    # If any combinations lack both FASTQ and BAM, report and exit
    if all_invalid_paths:
        print("No FASTQ or BAM found for the following combinations:")
        for read1, read2, bam in all_invalid_paths:
            print(f"  Read1: {read1}")
            print(f"  Read2: {read2}")
            print(f"  BAM  : {bam}")
            print("  ---")
        sys.exit(1)

    # Write CSVs with your mapping order (bam is in the mapping, right after seq_type)
    for sample_casename, rows in grouped_data.items():
        if sample_casename[1].startswith("patient_"):
            pass
        output_file = os.path.join(
            os.getcwd(), f"{sample_casename[0]}_{sample_casename[1]}.csv"
        )
        with open(output_file, "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=column_mapping.values())
            writer.writeheader()
            # ensure every row has bam field (and all others)
            for r in rows:
                for k in column_mapping.values():
                    r.setdefault(k, "")
                writer.writerow(r)
        print(output_file)


# column mapping (bam placed immediately after seq_type)
column_mapping = {
    "Patient ID": "sample",
    "Library ID": "library",
    "read1": "read1",
    "read2": "read2",
    "Enrichment step": "sample_captures",
    "Matched RNA-seq lib": "Matched_RNA",
    "Matched normal": "Matched_normal",
    "Diagnosis": "Diagnosis",
    "Case Name": "casename",
    "Type": "type",
    "FCID": "FCID",
    "Type of sequencing": "seq_type",
    "bam": "bam",
    "SampleRef": "genome",
}

# Call the function to process all samplesheets
process_samplesheets(samplesheet_dir, inputdir, column_mapping, sample_id, case_name)
