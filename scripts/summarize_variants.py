#!/usr/bin/env python3
"""
Summarise Mpox variant calling results from iVar TSV files.

This script walks through a directory that contains one subfolder per sample
(e.g. sample_S01, sample_S02, ...). Inside each sample_* folder it expects a
file called `variants.tsv` produced by the iVar variant-calling step.

For each sample, the script:
  - counts the total number of variant records,
  - counts how many variants have ALT_QUAL ≥ 20 and ≥ 30,
  - calculates the percentage of variants passing these quality thresholds,
  - calculates the transition/transversion (Ts/Tv) ratio for simple SNVs,
  - writes a tab-delimited summary table with one row per sample.

This summary will be used for downstream QC plots and to compare variant
quality across different pooling strategies or sequencing runs.
"""

import os
import argparse
import pandas as pd


def parse_variant_file(filepath):
    """
    Parse a single `variants.tsv` file and compute summary metrics.

    Parameters
    ----------
    filepath : str
        Path to a tab-delimited iVar variant file.

    Returns
    -------
    total : int
        Total number of variants in the file.
    qual20 : int
        Number of variants with ALT_QUAL ≥ 20.
    qual30 : int
        Number of variants with ALT_QUAL ≥ 30.
    pct20 : float
        Percentage of variants with ALT_QUAL ≥ 20.
    pct30 : float
        Percentage of variants with ALT_QUAL ≥ 30.
    tstv : float or str
        Transition/transversion ratio for simple SNVs (A,C,G,T only),
        or 'NA' if not calculable (e.g. no transversions).
    """
    try:
        # Read the TSV file produced by iVar.
        df = pd.read_csv(filepath, sep='\t')
    except Exception:
        # If file is missing, corrupted, or unreadable, return zeros/NA.
        return 0, 0, 0, 0.0, 0.0, 'NA'

    # If the file is empty, there are no variants to summarise.
    if df.empty or df.shape[0] == 0:
        return 0, 0, 0, 0.0, 0.0, 'NA'

    # Ensure that the expected columns are present.
    # ALT_QUAL: quality of the alternative allele;
    # REF/ALT: reference and alternate bases.
    required_cols = {'ALT_QUAL', 'REF', 'ALT'}
    if not required_cols.issubset(df.columns):
        return 0, 0, 0, 0.0, 0.0, 'NA'

    # Total number of variant records (all sites).
    total = df.shape[0]

    # Count how many variants pass quality thresholds (ALT_QUAL ≥ 20, ≥ 30).
    qual20 = df[df['ALT_QUAL'] >= 20].shape[0]
    qual30 = df[df['ALT_QUAL'] >= 30].shape[0]

    # Convert counts into percentages of all variants.
    pct20 = round((qual20 / total) * 100, 2) if total > 0 else 0.0
    pct30 = round((qual30 / total) * 100, 2) if total > 0 else 0.0

    # Restrict to simple single-nucleotide variants (SNVs) where both
    # REF and ALT are a single base and belong to {A, C, G, T}.
    snps = df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
    snps = snps[
        snps['REF'].isin(['A', 'T', 'C', 'G']) &
        snps['ALT'].isin(['A', 'T', 'C', 'G'])
    ]

    # Transitions: A<->G or C<->T changes.
    transitions = snps[
        ((snps['REF'] == 'A') & (snps['ALT'] == 'G')) |
        ((snps['REF'] == 'G') & (snps['ALT'] == 'A')) |
        ((snps['REF'] == 'C') & (snps['ALT'] == 'T')) |
        ((snps['REF'] == 'T') & (snps['ALT'] == 'C'))
    ].shape[0]

    # Transversions: all other SNV changes (more disruptive biologically).
    transversions = snps.shape[0] - transitions

    # Ts/Tv ratio: common QC metric for variant calls.
    # If there are no transversions, we return 'NA' to avoid division by zero.
    tstv = round(transitions / transversions, 2) if transversions > 0 else 'NA'

    return total, qual20, qual30, pct20, pct30, tstv


def main(input_dir, output_file):
    """
    Walk over all sample_* subdirectories, summarise their variants.tsv,
    and write a combined summary table.

    Parameters
    ----------
    input_dir : str
        Path to a directory containing subfolders named `sample_*`.
        Each subfolder is expected to contain a `variants.tsv` file.
    output_file : str
        Path to the output TSV file that will store the summary.
    """
    summary = []

    # Iterate over everything in the input directory.
    for entry in os.listdir(input_dir):
        full_path = os.path.join(input_dir, entry)

        # Only process subfolders whose name starts with "sample_".
        if os.path.isdir(full_path) and entry.startswith("sample_"):
            # Sample ID is taken from the folder name after "sample_".
            sample_id = entry.replace("sample_", "")

            # Expected variant file inside each sample folder.
            variant_file = os.path.join(full_path, "variants.tsv")

            if os.path.exists(variant_file):
                # Parse the variant file and compute metrics.
                total, q20, q30, pct20, pct30, tstv = parse_variant_file(
                    variant_file
                )
            else:
                # If no variants.tsv file exists, record zeros/NA for this sample.
                total, q20, q30, pct20, pct30, tstv = 0, 0, 0, 0.0, 0.0, 'NA'

            # Store one row per sample for the output table.
            summary.append([sample_id, total, q20, q30, pct20, pct30, tstv])

    # Convert the summary list into a DataFrame.
    df_out = pd.DataFrame(
        summary,
        columns=[
            "Sample",
            "TotalVariants",
            "VariantsQUAL20",
            "VariantsQUAL30",
            "PctQUAL20",
            "PctQUAL30",
            "TsTv"
        ]
    )

    # Write the combined summary as a tab-delimited file.
    df_out.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    # Command-line interface to run the script easily.
    parser = argparse.ArgumentParser(
        description=(
            "Summarise Mpox variant calling results from iVar TSV files. "
            "Expects an input folder containing sample_* subdirectories, "
            "each with a 'variants.tsv' file."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to folder with sample_* subfolders (each containing variants.tsv)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output summary TSV file."
    )
    args = parser.parse_args()

    main(args.input, args.output)
