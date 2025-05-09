import pandas as pd
import argparse

def main(args):
    data = pd.read_csv(args.proteomic_database,sep="\t")

    # ADD FILTERING LOGIC TO REMOVE CLUSTERS WITH TOO MANY VH GENES

    cluster_reps = {}
    for i in data['ClusterID'].unique():
        rep = data[data['ClusterID']==i]['CDRH3_aa'].value_counts().idxmax()
        cluster_reps[i] = rep
    with open(args.fasta,'w') as outf:
        for rep,sequence in cluster_reps.items():
            outf.write(f'>{rep}\n{sequence}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map Proteomic Database CSV to FASTA of cluster reps.')
    parser.add_argument('--proteomic_database', type=str, help='Proteomic Database CSV pointer')
    parser.add_argument('--fasta', type=str, help='Output FASTA filename')
    
    args = parser.parse_args()
    main(args)