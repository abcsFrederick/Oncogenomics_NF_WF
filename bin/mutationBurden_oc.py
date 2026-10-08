#!/usr/bin/env python3
"""Compute tumor mutation burden (nonsynonymous and all coding somatic calls).

Adapted from the COMPASS ngs_pipeline mutationBurden.py. Variant consequence is
classified using the OpenCRAVAT "Sequence Ontology" column instead of ANNOVAR's
ExonicFunc_refGene, since the annotatedFull files produced here are OpenCRAVAT-based.
"""
import sys

import pandas as pd

# Sequence Ontology terms counted as nonsynonymous (included in both nonsynonymous and all).
NONSYNONYMOUS_SO = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_elongation",
    "frameshift_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "complex_substitution",
    "splice_site_variant",
    "NMD_transcript_variant",
    "stop_retained_variant",
    "lnc_RNA",
}
# Silent coding term: included in "all somatic calls" but not in nonsynonymous.
SYNONYMOUS_SO = {"synonymous_variant"}
# Coding terms contributing to the "all somatic calls" count.
CODING_SO = NONSYNONYMOUS_SO | SYNONYMOUS_SO

# Normal block appears before the tumor block, so pandas suffixes the tumor
# duplicates with ".1".
NCOV_COL = "TotalCoverage"  # normal coverage
TCOV_COL = "TotalCoverage.1"  # tumor coverage
VAF_COL = "Variant Allele Freq.1"  # tumor VAF
SO_COL = "Sequence Ontology"


def get_total_bases(intervals_path):
    total_bp = 0
    with open(intervals_path, "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            total_bp += int(fields[2]) - int(fields[1])
    return total_bp


def main():
    if len(sys.argv) != 6:
        sys.exit(
            "Usage: mutationBurden_oc.py <ann_var> <intervals> <tcov> <ncov> <vaf>"
        )

    ann_var = sys.argv[1]
    intervals = sys.argv[2]
    tcov = float(sys.argv[3])
    ncov = float(sys.argv[4])
    vaf = float(sys.argv[5])

    total_bp = get_total_bases(intervals)
    if total_bp <= 0:
        sys.exit(f"No positive interval length found in: {intervals}")

    annvar = pd.read_csv(ann_var, sep="\t", dtype=str, keep_default_na=False)
    for col in (NCOV_COL, TCOV_COL, VAF_COL, SO_COL):
        if col not in annvar.columns:
            sys.exit(f"Missing required column '{col}' in {ann_var}")

    normal_cov = pd.to_numeric(annvar[NCOV_COL], errors="coerce")
    tumor_cov = pd.to_numeric(annvar[TCOV_COL], errors="coerce")
    tumor_vaf = pd.to_numeric(annvar[VAF_COL], errors="coerce")
    so_terms = annvar[SO_COL].str.strip()

    passes = (normal_cov >= ncov) & (tumor_cov >= tcov) & (tumor_vaf >= vaf)

    all_count = int((passes & so_terms.isin(CODING_SO)).sum())
    nonsyn_count = int((passes & so_terms.isin(NONSYNONYMOUS_SO)).sum())

    print(f"Total bases\t{total_bp}\n")
    print("Nonsynonymous somatic calls:")
    print(f"Mutation burden\t{nonsyn_count}")
    print(f"Mutation burden per megabase\t{(nonsyn_count / total_bp) * 1000000:.3f}\n")
    print("All somatic calls:")
    print(f"Mutation burden\t{all_count}")
    print(f"Mutation burden per megabase\t{(all_count / total_bp) * 1000000:.3f}\n")


if __name__ == "__main__":
    main()
