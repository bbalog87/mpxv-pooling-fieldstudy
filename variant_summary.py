#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def parse_variant_file(filepath):
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        return 0, 0, 0, 0.0, 0.0, 'NA'

    if df.empty or df.shape[0] == 0:
        return 0, 0, 0, 0.0, 0.0, 'NA'

    if 'ALT_QUAL' not in df.columns or 'REF' not in df.columns or 'ALT' not in df.columns:
        return 0, 0, 0, 0.0, 0.0, 'NA'

    total = df.shape[0]
    qual20 = df[df['ALT_QUAL'] >= 20].shape[0]
    qual30 = df[df['ALT_QUAL'] >= 30].shape[0]

    pct20 = round((qual20 / total) * 100, 2) if total > 0 else 0.0
    pct30 = round((qual30 / total) * 100, 2) if total > 0 else 0.0

    snps = df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
    snps = snps[snps['REF'].isin(['A', 'T', 'C', 'G']) & snps['ALT'].isin(['A', 'T', 'C', 'G'])]

    transitions = snps[((snps['REF'] == 'A') & (snps['ALT'] == 'G')) |
                       ((snps['REF'] == 'G') & (snps['ALT'] == 'A')) |
                       ((snps['REF'] == 'C') & (snps['ALT'] == 'T')) |
                       ((snps['REF'] == 'T') & (snps['ALT'] == 'C'))].shape[0]

    transversions = snps.shape[0] - transitions
    tstv = round(transitions / transversions, 2) if transversions > 0 else 'NA'

    return total, qual20, qual30, pct20, pct30, tstv

def main(input_dir, output_file):
    summary = []

    for entry in os.listdir(input_dir):
        full_path = os.path.join(input_dir, entry)
        if os.path.isdir(full_path) and entry.startswith("sample_"):
            sample_id = entry.replace("sample_", "")
            variant_file = os.path.join(full_path, "variants.tsv")
            if os.path.exists(variant_file):
                total, q20, q30, pct20, pct30, tstv = parse_variant_file(variant_file)
            else:
                total, q20, q30, pct20, pct30, tstv = 0, 0, 0, 0.0, 0.0, 'NA'

            summary.append([sample_id, total, q20, q30, pct20, pct30, tstv])

    df_out = pd.DataFrame(summary, columns=[
        "Sample", "TotalVariants", "VariantsQUAL20", "VariantsQUAL30",
        "PctQUAL20", "PctQUAL30", "TsTv"
    ])
    df_out.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Mpox variant summary from iVar TSVs.")
    parser.add_argument("-i", "--input", required=True, help="Path to folder with sample_* subfolders.")
    parser.add_argument("-o", "--output", required=True, help="Path to output summary TSV.")
    args = parser.parse_args()

    main(args.input, args.output)
