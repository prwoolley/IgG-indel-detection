import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import argparse
import os


def make_cluster_df(matches, peak_ratio):
    # Make cluster_df
    cluster_df = {'target_id':[],'std':[],'num_seqs':[],'tlen':[], 'num_peaks':[],'peak_lengths':[],'peak_heights':[]}
    for target_id in matches['target_id'].unique():
        subset = matches[matches['target_id']==target_id]
        cluster_df['tlen'].append(subset['tlen'].iloc[0])
        subset = subset['qlen']
        cluster_df['target_id'].append(target_id)
        cluster_df['std'].append(np.std(subset))
        cluster_df['num_seqs'].append(len(subset))
        height_threshold = peak_ratio*len(subset)
        try:
            start = min(subset)-5
            end = max(subset)+5
            counts, pos = np.histogram(subset, bins=np.linspace(start,end,end-start+1))
            peaks, _ = find_peaks(counts, height=height_threshold)
            peak_lengths = '_'.join([str(int(pos[peak])) for peak in peaks])
            peak_heights = '_'.join([str(int(counts[peak])) for peak in peaks])
            cluster_df['num_peaks'].append(len(peaks))
            cluster_df['peak_lengths'].append(peak_lengths)
            cluster_df['peak_heights'].append(peak_heights)
        except:
            cluster_df['num_peaks'].append(0)
            cluster_df['peak_lengths'].append('NA')
            cluster_df['peak_heights'].append('NA')
    cluster_df = pd.DataFrame(cluster_df)
    return cluster_df


def make_pd_cluster_df(proteomic_database, peak_ratio):
    # Load proteomic database and make pd_cluster_df the same way as cluster_df
    proteomic_database = pd.read_csv(proteomic_database,header=None,skiprows=1,names=['query_id','target_id','CDRH1_aa','CDRH2_aa','CDRH3_aa','VH','DH','JH','V_AvgSHM','Hiso','qseq'],sep="\t")
    proteomic_database['qlen'] = proteomic_database['qseq'].apply(lambda x: len(str(x)))
    pd_cluster_df = {'target_id':[],'std':[],'num_seqs':[],'num_peaks':[],'peak_lengths':[],'peak_heights':[]}
    for target_id in proteomic_database['target_id'].unique():
        subset = proteomic_database[proteomic_database['target_id']==target_id]
        subset = subset['qlen']
        pd_cluster_df['target_id'].append(target_id)
        pd_cluster_df['std'].append(np.std(subset))
        pd_cluster_df['num_seqs'].append(len(subset))
        height_threshold = peak_ratio*len(subset)
        try:
            start = min(subset)-5
            end = max(subset)+5
            counts, pos = np.histogram(subset, bins=np.linspace(start,end,end-start+1))
            peaks, _ = find_peaks(counts, height=height_threshold)
            peak_lengths = '_'.join([str(int(pos[peak])) for peak in peaks])
            peak_heights = '_'.join([str(int(counts[peak])) for peak in peaks])
            pd_cluster_df['num_peaks'].append(len(peaks))
            pd_cluster_df['peak_lengths'].append(peak_lengths)
            pd_cluster_df['peak_heights'].append(peak_heights)
        except:
            pd_cluster_df['num_peaks'].append(0)
            pd_cluster_df['peak_lengths'].append('NA')
            pd_cluster_df['peak_heights'].append('NA')
    pd_cluster_df = pd.DataFrame(pd_cluster_df)
    return pd_cluster_df


def merge_clusters(cluster_df, pd_cluster_df):
    # Write clusters present in proteomic database but not in cluster_df
    dropped_clusters = pd_cluster_df['target_id'][~pd_cluster_df['target_id'].isin(cluster_df['target_id'])].values
    dropped_clusters = np.sort(dropped_clusters)
    proteomic_database = proteomic_database[~proteomic_database['target_id'].isin(dropped_clusters)]
    pd_cluster_df = pd_cluster_df[~pd_cluster_df['target_id'].isin(dropped_clusters)]

    # Compute clusters that had a increase in number of peaks
    merged_df = pd.merge(cluster_df[['target_id','num_peaks']],pd_cluster_df[['target_id','num_peaks']],on='target_id',how='inner')
    merged_df['peak_difference'] = merged_df['num_peaks_x'] - merged_df['num_peaks_y']
    cluster_df = pd.merge(cluster_df,merged_df[['target_id','peak_difference']],on='target_id',how='inner')
    return cluster_df


def filter_assembled_annotation(assembled_annotation_in, matches, cluster_df, minimum_cluster_size):
    assembled_annotation_in = pd.read_csv(assembled_annotation_in,sep='\t')
    subset = cluster_df[(cluster_df['num_seqs']>minimum_cluster_size) & (cluster_df['num_peaks']>1)]
    targets = subset['target_id']
    subset = matches[matches['target_id'].isin(targets)]
    subset_df = subset.loc[subset.groupby('query_id')['bit_score'].idxmax()]
    aa_file_subset = assembled_annotation_in[assembled_annotation_in['readId'].isin(subset_df['query_id'].astype(int))]
    aa_file_subset = pd.merge(assembled_annotation_in, subset_df[['query_id','target_id']], left_on='readId', right_on = 'query_id')
    return aa_file_subset


def main(args):
    peak_ratio = args.peak_ratio
    proteomic_database = args.proteomic_database
    output_dir = args.output_dir
    matches = args.matches
    assembled_annotation_in = args.assembled_annotation_in
    minimum_cluster_size = args.minimum_cluster_size

    # Load matches file
    matches = pd.read_csv(matches,sep='\t',header=None,
                          names=['query_id','target_id','percent_identity','alignment_length','mismatches','gapopen','qstart','qend','tstart','tend','evalue','bit_score','qlen','tlen'])
    matches = matches.loc[matches.groupby('query_id')['bit_score'].idxmax()] # Get best match for each query_id

    # Make dataframes
    cluster_df = make_cluster_df(matches, peak_ratio)
    pd_cluster_df = make_pd_cluster_df(proteomic_database, peak_ratio)
    cluster_df = merge_clusters(cluster_df, pd_cluster_df)
    aa_file_subset = filter_assembled_annotation(assembled_annotation_in, matches, cluster_df, minimum_cluster_size)

    # Write output files
    cluster_df.to_csv(os.path.join(output_dir,'cluster_stats.csv'),index=False)
    aa_file_subset.to_csv(os.path.join(output_dir,f'{assembled_annotation_in.split('/')[-1]}_processed'),sep=',',index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Postprocess MMSeqs2 clusters. Outputs an text file of dropped clusters from original Proteomic Database, a csv of cluster statistics, and an abbreviated assembled_annotations text file.')
    parser.add_argument('--matches', type=str, help='MMSeqs2 matches file')
    parser.add_argument('--minimum_cluster_size', type=int, default=1000, help='Minimum size of a cluster to consider. Default is 1000')
    parser.add_argument('--peak_ratio', type=float, help='Ratio of sequences in a cluster that must be a peak. Default=0.1', default=0.1)
    parser.add_argument('--proteomic_database', type=str, help='Proteomic Database CSV pointer')
    parser.add_argument('--assembled_annotation_in', type=str, help='Assembled annotation file to subset')
    parser.add_argument('--output_dir', type=str, help='Output directory')
    args = parser.parse_args()
    main(args)

    
    