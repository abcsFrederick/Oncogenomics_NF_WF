#!/usr/bin/env python3
"""Combine mutect SNP and strelka indel tumor mutation burden (TMB).

Adapted from the COMPASS ngs_pipeline combineTMB.py. Reads two mutationburden.txt
files (produced by mutationBurden.py) and sums their nonsynonymous and all
somatic call counts, reporting the combined burden per megabase.
"""
import sys


def read_burden(path):
    """Return (total_bp, nonsynonymous, all) counts from a mutationburden.txt file."""
    with open(path) as fh:
        lines = fh.readlines()
    if len(lines) < 8:
        sys.exit(f"Unexpected format (fewer than 8 lines) in: {path}")
    total_bp = lines[0].split("\t")[1].strip()
    nonsyn = int(lines[3].split("\t")[1].strip())
    all_calls = int(lines[7].split("\t")[1].strip())
    return total_bp, nonsyn, all_calls


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: combineTMB_oc.py <mutect_mutationburden.txt> <strelka_indels_mutationburden.txt>"
        )

    mutect_file = sys.argv[1]
    strelka_file = sys.argv[2]

    total_bp, mut_non, mut_all = read_burden(mutect_file)
    _, strel_non, strel_all = read_burden(strelka_file)

    total_bp_val = float(total_bp)
    if total_bp_val <= 0:
        sys.exit(f"No positive total bases found in: {mutect_file}")

    non = mut_non + strel_non
    all_calls = mut_all + strel_all

    print(f"Total bases\t{total_bp}\n")
    print("Nonsynonymous somatic calls:")
    print(f"Mutation burden\t{non}")
    print(f"Mutation burden per megabase\t{(non / total_bp_val) * 1000000:.3f}\n")
    print("All somatic calls:")
    print(f"Mutation burden\t{all_calls}")
    print(f"Mutation burden per megabase\t{(all_calls / total_bp_val) * 1000000:.3f}\n")


if __name__ == "__main__":
    main()
