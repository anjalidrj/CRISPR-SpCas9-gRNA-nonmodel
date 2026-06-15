#!/usr/bin/env python3
"""
SpCas9 Guide RNA Design Pipeline
Uses GuideScan2 for off-target scoring and RS3 for on-target scoring
(Hsu2013 and Chen2013 tracrRNA versions)

Usage:
    python scripts/pipeline.py \
        --fasta data/example/pp.fasta \
        --gff3 data/example/pp.gff3 \
        --bed data/example/cymC.bed \
        --index data/example/pp_genome \
        --tracr both \
        --outdir results/
"""

import argparse
import os
import subprocess
import sys

import pandas as pd
from Bio import SeqIO

# ── Argument parser ───────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description='SpCas9 guide RNA design pipeline'
    )
    parser.add_argument('--fasta',      required=True, help='Genome FASTA file')
    parser.add_argument('--gff3',       required=True, help='Genome GFF3 file')
    parser.add_argument('--bed',        required=True, help='Target gene BED file')
    parser.add_argument('--index',      required=True, help='GuideScan2 index prefix')
    parser.add_argument('--tracr',      default='both',
                        choices=['Hsu2013', 'Chen2013', 'both'],
                        help='tracrRNA version for RS3 scoring')
    parser.add_argument('--outdir',     default='results/',
                        help='Output directory')
    parser.add_argument('--mismatches', default=3, type=int,
                        help='Max mismatches for GuideScan2 (default: 3)')
    parser.add_argument('--threshold',  default=1, type=int,
                        help='GuideScan2 off-target threshold (default: 1)')
    return parser.parse_args()

# ── Step 1: Extract guide sequences ──────────────────────────────────────────

def extract_guides(fasta_path, bed_path):
    """Extract all NGG protospacers from the target BED region."""
    print('\n[Step 1] Extracting guide sequences from target region...')

    with open(bed_path) as f:
        line = f.readline().strip().split('\t')
    chrom  = line[0]
    start  = int(line[1])
    end    = int(line[2])
    gene   = line[3] if len(line) > 3 else 'gene'
    strand = line[5] if len(line) > 5 else '+'

    print(f'  Target: {gene} | {chrom}:{start}-{end} ({strand})')

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))
    if chrom not in genome:
        sys.exit(f'ERROR: {chrom} not found in FASTA')

    region_seq = genome[chrom].seq[start:end]
    print(f'  Region length: {len(region_seq)} bp')

    guides = []

    # Forward strand NGG
    for i in range(len(region_seq) - 22):
        protospacer = region_seq[i:i+20]
        pam         = region_seq[i+20:i+23]
        if str(pam)[1:] == 'GG' and 'N' not in str(protospacer).upper():
            guides.append({
                'guide_sequence' : str(protospacer).upper(),
                'chromosome'     : chrom,
                'start'          : start + i,
                'end'            : start + i + 20,
                'strand'         : '+',
                'pam'            : str(pam).upper(),
                'gene'           : gene
            })

    # Reverse strand NGG
    rev_seq = region_seq.reverse_complement()
    rev_len = len(rev_seq)
    for i in range(rev_len - 22):
        protospacer = rev_seq[i:i+20]
        pam         = rev_seq[i+20:i+23]
        if str(pam)[1:] == 'GG' and 'N' not in str(protospacer).upper():
            guides.append({
                'guide_sequence' : str(protospacer).upper(),
                'chromosome'     : chrom,
                'start'          : end - (i + 20),
                'end'            : end - i,
                'strand'         : '-',
                'pam'            : str(pam).upper(),
                'gene'           : gene
            })

    df = pd.DataFrame(guides)
    print(f'  Found {len(df)} candidate guides '
          f'({len(df[df.strand=="+"])} forward, '
          f'{len(df[df.strand=="-"])} reverse)')
    return df

# ── Step 2: Write kmers file ──────────────────────────────────────────────────

def write_kmers_file(guides_df, outdir):
    """Write guide sequences to a CSV file for GuideScan2 enumerate."""
    kmers_path = os.path.join(outdir, 'kmers.txt')
    with open(kmers_path, 'w') as f:
        f.write('id,sequence,pam,chromosome,position,sense\n')
        for i, row in guides_df.iterrows():
            f.write(f'{i},{row["guide_sequence"]},{row["pam"]},'
                    f'{row["chromosome"]},{row["start"]},{row["strand"]}\n')
    print(f'  Written {len(guides_df)} kmers to {kmers_path}')
    return kmers_path

# ── Step 3: Run GuideScan2 ────────────────────────────────────────────────────

def run_guidescan(index_prefix, kmers_path, outdir, mismatches=3, threshold=1):
    """Run GuideScan2 enumerate to score guides for off-target activity."""
    print('\n[Step 2] Running GuideScan2 off-target scoring...')

    gs2_output = os.path.join(outdir, 'guidescan2_output.csv')

    cmd = [
        'guidescan', 'enumerate',
        '--mismatches',  str(mismatches),
        '--threshold',   str(threshold),
        '--format',      'csv',
        '--mode',        'succinct',
        '--kmers-file',  kmers_path,
        '--output',      gs2_output,
        index_prefix
    ]

    print(f'  Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print('ERROR running GuideScan2:')
        print(result.stderr)
        sys.exit(1)

    print(f'  GuideScan2 complete. Output: {gs2_output}')
    return gs2_output

# ── Step 4: Parse GuideScan2 output ──────────────────────────────────────────

def parse_guidescan_output(gs2_output, guides_df):
    """Parse GuideScan2 CSV and merge specificity scores with guides."""
    print('\n[Step 3] Parsing GuideScan2 output...')

    gs2_df = pd.read_csv(gs2_output)

    # Keep only exact matches (distance=0) to get one row per guide
    gs2_exact = gs2_df[gs2_df['match_distance'] == 0].copy()
    gs2_exact = gs2_exact[['id', 'specificity']].drop_duplicates('id')

    # Merge with guides
    guides_df = guides_df.reset_index(drop=True)
    guides_df['id'] = guides_df.index
    merged = guides_df.merge(gs2_exact, on='id', how='left')
    merged['specificity'] = merged['specificity'].fillna(0)

    print(f'  Merged {len(merged)} guides with GuideScan2 scores')
    print(f'  Specificity range: '
          f'{merged["specificity"].min():.3f} – '
          f'{merged["specificity"].max():.3f}')
    return merged

# ── Step 5: RS3 on-target scoring ─────────────────────────────────────────────

def score_guides_rs3(merged_df, fasta_path, tracr='both'):
    """Score guides with RS3 (Hsu2013 and/or Chen2013 tracrRNA).
    RS3 requires 30nt context: 4nt upstream + 20nt guide + 3nt PAM + 3nt downstream.
    """
    print('\n[Step 4] Running RS3 on-target scoring...')

    from rs3.seq import predict_seq
    import numpy as np

    # Load genome to get flanking sequences
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))

    sequences_30nt = []
    for _, row in merged_df.iterrows():
        chrom  = row['chromosome']
        start  = int(row['start'])
        end    = int(row['end'])
        strand = row['strand']
        seq    = genome[chrom].seq

        if strand == '+':
            # 4nt upstream + 20nt guide + 3nt PAM + 3nt downstream
            context = seq[start-4 : end+6]
            context_str = str(context).upper()
        else:
            # reverse complement: 3nt upstream of PAM + 3nt PAM + 20nt guide + 4nt downstream
            context = seq[start-6 : end+4]
            context_str = str(context.reverse_complement()).upper()

        if len(context_str) == 30:
            sequences_30nt.append(context_str)
        else:
            sequences_30nt.append('N' * 30)  # edge case fallback

    tracr_versions = []
    if tracr in ['Hsu2013', 'both']:
        tracr_versions.append('Hsu2013')
    if tracr in ['Chen2013', 'both']:
        tracr_versions.append('Chen2013')

    for tracr_v in tracr_versions:
        print(f'  Scoring with tracrRNA: {tracr_v}...')
        scores = predict_seq(sequences_30nt, sequence_tracr=tracr_v)
        scores = np.array(scores)
        merged_df[f'rs3_{tracr_v}'] = scores
        print(f'  Score range: {scores.min():.3f} – {scores.max():.3f}')

    return merged_df

# ── Step 6: Rank guides ───────────────────────────────────────────────────────

def rank_guides(scored_df, tracr='both'):
    """Rank guides by specificity then RS3 score."""
    print('\n[Step 5] Ranking guides...')

    # Sort by specificity (desc) then RS3 Chen2013 or Hsu2013 (desc)
    sort_cols = ['specificity']
    if 'rs3_Chen2013' in scored_df.columns:
        sort_cols.append('rs3_Chen2013')
    elif 'rs3_Hsu2013' in scored_df.columns:
        sort_cols.append('rs3_Hsu2013')

    ranked_df = scored_df.sort_values(sort_cols, ascending=False)
    ranked_df = ranked_df.reset_index(drop=True)
    ranked_df.index += 1
    ranked_df.index.name = 'rank'

    print(f'  Top 5 guides:')
    display_cols = ['guide_sequence', 'chromosome', 'start', 'end',
                    'strand', 'specificity'] + \
                   [c for c in ranked_df.columns if c.startswith('rs3_')]
    print(ranked_df[display_cols].head().to_string())
    return ranked_df

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Step 1 — Extract guides
    guides_df  = extract_guides(args.fasta, args.bed)
    raw_path   = os.path.join(args.outdir, 'guides_raw.csv')
    guides_df.to_csv(raw_path, index=False)
    print(f'  Saved: {raw_path}')

    # Step 2 — Write kmers file
    kmers_path = write_kmers_file(guides_df, args.outdir)

    # Step 3 — GuideScan2 off-target scoring
    gs2_output = run_guidescan(
        args.index, kmers_path, args.outdir,
        args.mismatches, args.threshold
    )

    # Step 4 — Parse GuideScan2 output and merge
    merged_df  = parse_guidescan_output(gs2_output, guides_df)

    # Step 5 — RS3 on-target scoring
    scored_df  = score_guides_rs3(merged_df, args.fasta, args.tracr)
    scored_path = os.path.join(args.outdir, 'guides_scored.csv')
    scored_df.to_csv(scored_path, index=False)
    print(f'  Saved: {scored_path}')

    # Step 6 — Rank and save final output
    ranked_df  = rank_guides(scored_df, args.tracr)
    ranked_path = os.path.join(args.outdir, 'guides_ranked.csv')
    ranked_df.to_csv(ranked_path)
    print(f'  Saved: {ranked_path}')

    print('\n✅ Pipeline complete!')
    print(f'   Raw guides:    {raw_path}')
    print(f'   Scored guides: {scored_path}')
    print(f'   Ranked guides: {ranked_path}')