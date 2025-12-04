#!/usr/bin/env python3

"""
Usage:
  python3 maf_to_alt_ref.py --maf all.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_allColumns.txt \
                           --fasta reference.fa \
                           --out_prefix output_prefix

This script:
 - Selects columns (1-based) in this order: c12, c3, c15, c20, c19, c13, c14, c31
 - Renames them to: Chr, Pos, Type, Coverage, Normal_MAF, Tumor_MAF, Ref_Allele, Alt_Allele, Score
 - Inserts a new column "Size" at position 4 (1-based) with all values = 1
 - Inserts a new column "Chr" at position 1 (1-based) with all values = 'chrM'
 - Writes a trimmed TSV and creates an alternate FASTA with all variants applied (non-overlapping required).
"""

import argparse
import sys
import pandas as pd
from collections import defaultdict

try:
    from Bio import SeqIO
except Exception:
    SeqIO = None

def read_maf(maf_path, selected_cols):
    # Read with header if present; use pandas and then select by index (1-based -> 0-based)
    df = pd.read_csv(maf_path, sep='\t', low_memory=False, index_col=False)
    # If file has fewer columns than requested indices, raise
    max_idx = max(selected_cols)
    if df.shape[1] < max_idx:
        raise SystemExit(f"MAF file has {df.shape[1]} columns but requested column index {max_idx}.")
    # select by position
    selected = [df.iloc[:, i-1].copy() for i in selected_cols]
    return pd.concat(selected, axis=1)

def make_alt_ref(fasta_path, variants_df, out_fasta_path):
    # variants_df columns: Chr, Pos, Ref_Allele, Alt_Allele
    # Group by chromosome
    chr_groups = defaultdict(list)
    for _, r in variants_df.iterrows():
        chr_name = str('chrM')
        pos = int(r['Pos'])
        ref = str(r['Ref_Allele'])
        alt = str(r['Alt_Allele'])
        start = pos - 1 # convert to 0-based for string slicing later
        end = start + len(ref)
        chr_groups[chr_name].append({'start': start, 'end': end, 'ref': ref, 'alt': alt, 'pos': pos})

    # Read fasta sequences
    seqs = {}
    if SeqIO:
        for rec in SeqIO.parse(fasta_path, 'fasta'):
            seqs[rec.id] = str(rec.seq)
    else:
        # simple fasta parser
        with open(fasta_path) as fh:
            cur_id = None
            cur_seq = []
            for line in fh:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if cur_id:
                        seqs[cur_id] = ''.join(cur_seq)
                    cur_id = line[1:].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_id:
                seqs[cur_id] = ''.join(cur_seq)

    # Convert seq chromosome 'MT' to 'chrM' if needed
    if 'MT' in seqs and 'chrM' not in seqs:
        seqs['chrM'] = seqs.pop('MT')
    
    # Build alternate sequences
    out_records = []
    for chrom, seq in seqs.items():
        chrom_variants = chr_groups.get(chrom, [])
        if not chrom_variants:
            out_records.append((chrom, seq))
            continue

        # validate and sort
        chrom_variants.sort(key=lambda x: x['start'])
        # check overlaps
        prev_end = -1
        for v in chrom_variants:
            if v['start'] < prev_end:
                print(f"Overlapping variants detected on {chrom} at pos {v['pos']}. Skipping to next variant.")
                continue
            # validate reference matches
            ref_segment = seq[v['start']:v['end']]
            if ref_segment.upper() != v['ref'].upper():
                # allow replacement if ref allele doesn't match: warn and continue (use provided ref)
                print(f"Warning: reference mismatch on {chrom}:{v['pos']} (expected '{v['ref']}' vs fasta '{ref_segment}'). Using provided ref.", file=sys.stderr)
            prev_end = v['end']

        # apply variants
        out_seq_parts = []
        cursor = 0
        prev_end = -1
        for v in chrom_variants:
            if v['start'] < prev_end:
                print(f"Overlapping variants detected on {chrom} at pos {v['pos']}. Skipping to next variant.")
                continue
            out_seq_parts.append(seq[cursor:v['start']])
            out_seq_parts.append(v['alt'])
            cursor = v['end']
        out_seq_parts.append(seq[cursor:])
        new_seq = ''.join(out_seq_parts)
        out_records.append((chrom, new_seq))

    # # Do pairwise alignment of the original fasta sequence and new sequence to verify correctness
    # from Bio import Align
    # aligner = Align.PairwiseAligner()
    # aligner.mode = 'global'
    # aligner.open_gap_score = -10.0
    # aligner.extend_gap_score = -5.0
    # aligner.match_score = 1.0
    # aligner.mismatch_score = -1.0
    # aligner.end_extend_gap_score = -0.5
    # alignments = aligner.align(seq, out_records[0][1])
    # top_align = alignments[0]
    # with open(f"text_alignment.txt", 'w') as align_fh:
    #     print(top_align, file = align_fh)

    # write out FASTA
    with open(out_fasta_path, 'w') as outfh:
        for chrom, seq in out_records:
            outfh.write(f">{chrom}\n")
            # wrap at 60 chars per line
            for i in range(0, len(seq), 70):
                outfh.write(seq[i:i+70] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Trim MAF-like file to selected columns and create alternate FASTA with variants applied.")
    parser.add_argument('--maf', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--out_prefix', required=True)
    args = parser.parse_args()

    # Select columns (1-based): c13, c4, c16, c21, c20, c14, c15, c32
    selected_cols = [12, 3, 15, 20, 19, 13, 14, 31]  # 0-based indices
    new_names     = ['Pos', 'Type', 'Coverage', 'Normal_MAF', 'Tumor_MAF', 'Ref_Allele', 'Alt_Allele', 'Score']

    if len(selected_cols) != len(new_names):
        raise SystemExit("selected_cols and new_names length mismatch. Edit the lists in the script.")

    # read and select
    sel_df = read_maf(args.maf, selected_cols)
    sel_df.columns = new_names

    # Insert Chr column at position 1 (1-based) => index 0 (0-based), fill with 'chrM'
    sel_df.insert(0, 'Chr', 'chrM')

    # Insert Size column at position 4 (1-based) => index 3 (0-based), fill with '1'
    sel_df.insert(3, 'Size', '1')

    # ensure Pos is integer
    sel_df['Pos'] = sel_df['Pos'].astype(int)

    # write trimmed TSV
    out_tsv = f"{args.out_prefix}_trimmed.tab"
    sel_df.to_csv(out_tsv, sep='\t', index=False)
    print(f"Wrote trimmed table to {out_tsv}")

    # Create alternate FASTA by applying all variants (non-overlapping) to reference
    # For applying, need columns Chr, Pos, Ref_Allele, Alt_Allele
    required = {'Chr', 'Pos', 'Ref_Allele', 'Alt_Allele'}
    if not required.issubset(set(sel_df.columns)):
        raise SystemExit(f"Need columns {required} to create alternate FASTA. Current columns: {list(sel_df.columns)}")

    out_fasta = f"{args.out_prefix}_alt_genome.fa"
    make_alt_ref(args.fasta, sel_df[['Chr','Pos','Ref_Allele','Alt_Allele']], out_fasta)
    print(f"Wrote alternate FASTA to {out_fasta}")

if __name__ == '__main__':
    main()