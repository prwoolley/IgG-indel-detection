# IgG-indel-detection
Indel detection software designed to be used on MiXCR outputs

Inputs are two files: Assembled-annotation and Proteomic-database

The pipeline is outlined as follows:

0. Creates directories output/fastas, output/mmseqs, output/igblast, where "output" is a user defined directory.
1. Pull cluster rep sequences from Proteomic-database and write to FASTA file.
2. Pull sequences from Assembled-annotation and write to FASTA file.
3. Create MMseqs database for the Assembled-annotation FASTA.
4. Create MMseqs database for the Proteomic-database FASTA.
5. Search Assembled-annotation database against Proteomic-database. Write alignments to m8 file.
6. Load m8 file and original Proteomic-databases. For each cluster, compute the number of peaks along the sequence length axis. Flag any clusters of a certain cluster size threshold that had an increase in the number of peaks. Open the Assembled-annotation file. Write any sequence from the flagged clusters to an abbreviated Assembled-annotation file and FASTA file.
7. Run IgBLAST on the abbreviated Assembled-annotation FASTA file.
8. Run insertion_finder.py on IgBLAST MSAs.