import pandas as pd
import argparse
from Bio.Seq import Seq


def translate_dna_to_protein(nucleotide_sequence):
    try:
        dna_seq = Seq(nucleotide_sequence)
        protein_sequence = dna_seq.translate()
        if "*" in str(protein_sequence):
            return None
        return str(protein_sequence)
    except Exception as e:
        return f"Error: {e}"


def main(args):
    assembled_annotation = pd.read_csv(args.assembled_annotation_in,sep='\t')
    assembled_annotation['translation'] = assembled_annotation['readSequence'].apply(lambda x: translate_dna_to_protein(x))
    assembled_annotation.dropna(subset=['translation'], inplace=True)

    with open(args.fasta,'w') as outf:
        for i,seq in enumerate(assembled_annotation['translation']):
            outf.write(f'>{i}\n{seq}\n')

    assembled_annotation.to_csv(args.assembled_annotation_out,index=False,sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map Proteomic Database CSV to FASTA of cluster reps.')
    parser.add_argument('--assembled_annotation_in', type=str, help='Input assembled annotation file')
    parser.add_argument('--output_dir', type=str, help='Proteomic Database CSV pointer')
    parser.add_argument('--fasta', type=str, help='Output FASTA filename')
    parser.add_argument('--assembled_annotation_out', type=str, help='Output assembled annotation file')
    
    args = parser.parse_args()
    main(args)