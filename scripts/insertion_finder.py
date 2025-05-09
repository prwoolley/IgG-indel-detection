import argparse
import subprocess
import os
import pandas as pd
import re


def simple_process(v_call):
    # Split on comma and take the first instance
    first_call = v_call.split(',')[0]
    # Split on asterisk and take the first instance
    base_gene = first_call.split('*')[0]
    return base_gene


def find_gap(seq):
    # Look for a continuous gap of '-' of length greater than 3 (i.e. at least 4)
    match = re.search(r'(-{4,})', seq)
    if match:
        gap_length = len(match.group(0))
        # Convert positions to 1-indexed
        gap_start = match.start() + 1  
        gap_end = match.end()         # .end() already gives the last position in 1-indexed terms when counting the gap
        return pd.Series([gap_length, gap_start, gap_end])
    else:
        # If no gap is found, return None for each value
        return pd.Series([None, None, None])


def main(args):
    df_igb = pd.read_csv(base_filename + '_igblast.tsv', sep='\t')
    
    df_igb = df_igb[['sequence_id', 'sequence', 'sequence_aa', 'v_call', 'v_germline_alignment',
           'v_germline_alignment_aa','fwr1', 'fwr1_aa', 'cdr1', 'cdr1_aa', 
           'fwr2', 'fwr2_aa', 'cdr2','cdr2_aa', 'fwr3', 'fwr3_aa']]
   
    df_igb['VH'] = df_igb['v_call'].apply(simple_process)
    df_igb[['gap_length', 'gap_start', 'gap_end']] = df_igb['v_germline_alignment_aa'].apply(find_gap)
    
    df_igb.to_csv(base_filename+'_aligned.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Script to find insertions and output their positions along the VH gene.",
        epilog="""Example usage: python3 insertion_finder.py cluster_stats.csv matches.m8 vh.assembled_annotation"""
    )
    # Define a positional argument for the input file
    parser.add_argument('cluster_stats', help='Path to the cluster_stats.csv file')
    parser.add_argument('matches', help='Path to the matches.m8 file')
    parser.add_argument('vh_assembled_annotation', help='Path to the assembled annotation file')
    args = parser.parse_args()
    main(args)