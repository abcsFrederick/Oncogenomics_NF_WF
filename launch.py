#!/usr/bin/env python3
"""
Launcher for the Oncogenomics Nextflow workflow (Biowulf/SLURM).

Usage (existing):
  launch.py --samplesheet /path/to/sheet.csv [--outdir DIR] [--genome {hg19,mm39}] [--platform biowulf]
            [--profile PROFILE] [--no-resume] [--no-cleanup]

Usage (new, auto-generate samplesheet first):
  launch.py --patient P123 --casename CASE_X [--inputdir /data/khanlab/projects/DATA] \
            [--outdir DIR] [--genome {hg19,hg38,mm39}] [--platform biowulf] [--profile PROFILE] \
            [--no-resume] [--no-cleanup]
"""
import argparse
import csv
import sys
import shlex
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Optional

# ---- Paths & defaults ----
DEFAULT_WF_HOME = Path("/data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF")
DEV_WF_HOME = Path("/data/khanlab/projects/Nextflow_dev/dev/vg_dev/Oncogenomics_NF_WF")

# CONFIG_FILE = WF_HOME / "nextflow.config"
DEFAULT_OUTDIR = "/data/khanlab/projects/processed_DATA"
DEFAULT_GENOME = "hg19"
DEFAULT_PLATFORM = "biowulf"
ACCEPTED_GENOMES = {"hg19", "hg38", "mm39"}

# Fixed builder script; users do NOT override this on CLI.
# Invoked as: python samplesheet_builder.py Patient Casename [inputdir]
# BUILDER_SCRIPT = WF_HOME / "bin" / "samplesheet_builder.py"


# ----------------- helpers -----------------
def read_unique_values(csv_path: Path, field: str) -> List[str]:
    vals: List[str] = []
    with csv_path.open(newline="") as fh:
        r = csv.DictReader(fh)
        if not r.fieldnames or field not in r.fieldnames:
            sys.exit(f"ERROR: required column '{field}' not found in {csv_path}")
        for row in r:
            v = (row.get(field) or "").strip()
            if v:
                vals.append(v)
    return sorted(set(vals))


def collapse_unique(vals: List[str], label: str) -> str:
    if not vals:
        sys.exit(f"ERROR: No non-empty values found for '{label}' in samplesheet.")
    if len(vals) > 1:
        sys.exit(
            f"ERROR: Multiple unique values found for '{label}' in samplesheet: {', '.join(vals)}"
        )
    return vals[0]


def run_cmd_argv(argv: List[str], cwd: Optional[Path] = None) -> None:
    """Run a command (argv form) with clear errors. If cwd is given, run in that directory."""
    try:
        proc = subprocess.run(
            argv,
            check=True,
            text=True,
            capture_output=True,
            cwd=str(cwd) if cwd else None,
        )
        if proc.stdout:
            print(proc.stdout.strip())
        if proc.stderr:
            # Show builder warnings/info but don't fail (return code already checked)
            sys.stderr.write(proc.stderr)
    except subprocess.CalledProcessError as e:
        msg = []
        msg.append("ERROR: command failed.")
        msg.append(f"Command: {' '.join(map(shlex.quote, argv))}")
        if e.stdout:
            msg.append("\n--- STDOUT ---\n" + e.stdout)
        if e.stderr:
            msg.append("\n--- STDERR ---\n" + e.stderr)
        sys.exit("\n".join(msg))


def find_samplesheet_after_build(resultsdir: Path, patient: str, casename: str) -> Path:
    """
    Expected: resultsdir/<patient>_<casename>.csv
    Fallback: if exactly one CSV exists directly under resultsdir, use it.
    """
    expected = resultsdir / f"{patient}_{casename}.csv"
    if expected.exists():
        return expected

    # Fallback: look for a single CSV file under resultsdir (non-recursive)
    csvs = sorted(resultsdir.glob("*.csv"))
    if len(csvs) == 1:
        print(f"ℹ️  Using detected samplesheet: {csvs[0]}")
        return csvs[0]
    elif len(csvs) == 0:
        sys.exit(
            "ERROR: Expected samplesheet not found after builder.\n"
            f"Tried: {expected}\n"
            f"No CSVs found in: {resultsdir}\n"
            "Please ensure samplesheet_builder.py writes the CSV with the expected name, "
            "or run with --samplesheet."
        )
    else:
        listing = "\n".join(str(p) for p in csvs)
        sys.exit(
            "ERROR: Multiple CSVs found; cannot infer the samplesheet uniquely.\n"
            f"Directory: {resultsdir}\nFound:\n{listing}\n"
            "Please remove extras or pass --samplesheet explicitly."
        )


# ----------------- CLI -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Submit the Oncogenomics Nextflow workflow to SLURM."
    )

    # Mutually exclusive source of samplesheet
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--samplesheet", help="Path to existing samplesheet CSV")
    src.add_argument(
        "--patient", help="Patient ID (requires --casename to build a samplesheet)"
    )

    # If --patient is used, require --casename too
    p.add_argument("--casename", help="Case name (required with --patient)")

    # Optional input directory for the builder:
    # The builder will be invoked as: python samplesheet_builder.py Patient Casename [inputdir]
    p.add_argument(
        "--inputdir",
        default=None,
        help="Optional input directory passed to samplesheet_builder.py",
    )

    # Existing launcher flags
    p.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=f"Results root directory (default: {DEFAULT_OUTDIR})",
    )
    p.add_argument(
        "--genome",
        default=DEFAULT_GENOME,
        choices=sorted(ACCEPTED_GENOMES),
        help=f"Genome (default: {DEFAULT_GENOME})",
    )
    p.add_argument(
        "--platform",
        default=DEFAULT_PLATFORM,
        help=f"Platform label (default: {DEFAULT_PLATFORM})",
    )
    p.add_argument(
        "--profile",
        default=None,
        help="Explicit Nextflow profile (overrides genome→profile mapping)",
    )
    p.add_argument(
        "--no-resume", action="store_true", help="Run Nextflow without -resume"
    )

    # Cleanup: enabled by default. Use --no-cleanup to opt out.
    p.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Disable cleanup (default: cleanup runs after successful job)",
    )
    p.add_argument(
        "--wf-home",
        default=None,
        help="Override workflow home (path containing main.nf and nextflow.config)",
    )
    p.add_argument(
        "--dev",
        action="store_true",
        help=f"Use dev WF at {str(DEV_WF_HOME)}",
    )

    return p.parse_args()


# ----------------- main -----------------
def main():
    args = parse_args()
    outdir_explicit = any(
        arg == "--outdir" or arg.startswith("--outdir=") for arg in sys.argv[1:]
    )

    if args.dev:
        WF_HOME = DEV_WF_HOME
        DEFAULT_OUTDIR_DEV = "/data/khanlab/projects/Nextflow_dev/dev/vg_dev"
        if not outdir_explicit:
            args.outdir = DEFAULT_OUTDIR_DEV
        print(f"🔧 Dev mode enabled — using WF_HOME={WF_HOME} and OUTDIR={args.outdir}")
    else:
        WF_HOME = DEFAULT_WF_HOME

    CONFIG_FILE = WF_HOME / "nextflow.config"
    BUILDER_SCRIPT = WF_HOME / "bin" / "samplesheet_builder.py"
    # Validate mode
    if args.patient and not args.casename:
        sys.exit("ERROR: --patient requires --casename.")
    if args.casename and not args.patient:
        sys.exit("ERROR: --casename requires --patient.")

    outdir_root = Path(args.outdir).resolve()
    genome = args.genome
    platform = args.platform
    do_cleanup = not args.no_cleanup  # cleanup is ON by default

    if genome not in ACCEPTED_GENOMES:
        sys.exit(
            f"Invalid genome specified. Accepted values: {', '.join(sorted(ACCEPTED_GENOMES))}"
        )

    # Resolve samplesheet
    if args.samplesheet:
        samplesheet = Path(args.samplesheet).resolve()
        if not samplesheet.exists():
            sys.exit(f"ERROR: samplesheet not found: {samplesheet}")
        patient = collapse_unique(read_unique_values(samplesheet, "sample"), "sample")
        casename = collapse_unique(
            read_unique_values(samplesheet, "casename"), "casename"
        )
    else:
        # Build samplesheet first using fixed builder
        patient, casename = args.patient, args.casename

        # Ensure resultsdir exists before running the builder there
        resultsdir = (outdir_root / patient / casename).resolve()
        resultsdir.mkdir(parents=True, exist_ok=True)

        if not BUILDER_SCRIPT.exists():
            sys.exit(f"ERROR: builder not found at {BUILDER_SCRIPT}")

        # Call builder: python samplesheet_builder.py Patient Casename [inputdir]  (run inside resultsdir)
        builder_cmd = ["python", str(BUILDER_SCRIPT), patient, casename]
        if args.inputdir:
            builder_cmd.append(args.inputdir)
        builder_cmd.extend(["--genome", shlex.quote(genome)])

        print(
            "🔧 Generating samplesheet via:\n  "
            + " ".join(map(shlex.quote, builder_cmd))
        )
        run_cmd_argv(builder_cmd, cwd=resultsdir)

        # Locate resulting samplesheet and validate
        samplesheet = find_samplesheet_after_build(resultsdir, patient, casename)
        _ = collapse_unique(read_unique_values(samplesheet, "sample"), "sample")
        _ = collapse_unique(read_unique_values(samplesheet, "casename"), "casename")

    # Common dirs (now that we know patient/casename)
    resultsdir = (outdir_root / patient / casename).resolve()
    logdir = resultsdir / "log"
    nxf_home = resultsdir / ".nextflow"
    for pth in (resultsdir, logdir, nxf_home):
        pth.mkdir(parents=True, exist_ok=True)

    # Profile selection
    if args.profile:
        profile = args.profile
    else:
        profile = (
            "biowulf_test_run_slurm"
            if genome in ("hg19", "hg38")
            else "biowulf_mouse_RNA_slurm"
        )

    logname = Path(samplesheet).stem
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_name = f"run_{logname}_{timestamp}"
    sb_log_out = resultsdir / f"{logname}_%A_{timestamp}.out"
    job_log_expr = f"{resultsdir}/{logname}_${{SLURM_JOB_ID}}_{timestamp}.out"
    resume_flag = "" if args.no_resume else "-resume"

    # SBATCH script
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name="{logname}"
#SBATCH --output="{sb_log_out}"
#SBATCH --cpus-per-task=2
#SBATCH --mem=5g
#SBATCH --time=08-00:00:00

set -euo pipefail
cd "{resultsdir}"

logfile="{job_log_expr}"
export NXF_HOME="{nxf_home}"

module load nextflow/23.10.0 singularity graphviz

nextflow run -c {shlex.quote(str(CONFIG_FILE))} -profile {shlex.quote(profile)} \\
  --logdir {shlex.quote(str(logdir))} {shlex.quote(str(WF_HOME / 'main.nf'))} {resume_flag} \\
  --samplesheet {shlex.quote(str(samplesheet))} \\
  --resultsdir {shlex.quote(str(outdir_root))} \\
  --genome_v {shlex.quote(genome)} \\
  --platform {shlex.quote(platform)} \\
  -name {shlex.quote(run_name)}
"""

    if do_cleanup:
        sbatch_script += f"""
if [[ -f "$logfile" ]]; then
  if tail -n 200 "$logfile" | grep -q "Completed at:" && \\
     tail -n 200 "$logfile" | grep -q "Succeeded"; then
    echo "✅ Run succeeded — cleaning up work dirs for {run_name}"
    nextflow clean -f {shlex.quote(run_name)}
    rm -rf "{nxf_home}"
  else
    echo "⚠️  Run did not finish — skipping cleanup"
  fi
fi
"""

    # Submit
    try:
        proc = subprocess.run(
            ["sbatch"], input=sbatch_script, text=True, check=True, capture_output=True
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f"\nERROR: sbatch submission failed.\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}\n"
        )
        sys.exit(2)

    print(proc.stdout.strip())
    print(f"Patient:      {patient}")
    print(f"Casename:     {casename}")
    print(f"Samplesheet:  {samplesheet}")
    print(f"Results:      {resultsdir}")
    print(f"Run name:     {run_name}")
    print(f"Profile:      {profile}")
    print(f"Resume:       {'off' if args.no_resume else 'on'}")
    print(f"Cleanup:      {'enabled' if do_cleanup else 'disabled'}")


if __name__ == "__main__":
    main()
