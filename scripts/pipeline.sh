#!/bin/bash

# Default values
OUTPUT_DIR="./output"
PEAK_RATIO=0.1
THREADS=8
SCRIPTS_DIR="/stor/work/Ellington/bioml/Dalton_task/scripts"
MINCLUSTERSIZE=1000


# Parse named arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --scripts_dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        --proteomic_database)
            PROTEOMIC_DATABASE="$2"
            shift 2
            ;;
        --assembled_annotation)
            ASSEMBLED_ANNOTATION="$2"
            shift 2
            ;;
        --igblast_path)
            IGBLAST_PATH="$2"
            shift 2
            ;;
        --peak_ratio)
            PEAK_RATIO="$2"
            shift 2
            ;;
        --minimum_cluster_size)
            MINCLUSTERSIZE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: pipeline.sh --output_dir <directory> --proteomic_database <file> --assembled_annotation <file> --scripts_dir <directory> [--peak_ratio <0.1 default>] [--threads <8 default>]"
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$PROTEOMIC_DATABASE" ] || [ -z "$ASSEMBLED_ANNOTATION" ] || [ -z "$IGBLAST_PATH"]; then
    echo "Error: --proteomic_database, --assembled_annotation, and --igblast_path parameters are required"
    exit 1
fi


# Create output directories if they do not exist
[ -d "${OUTPUT_DIR}" ] || mkdir "${OUTPUT_DIR}"
[ -d "${OUTPUT_DIR}/fastas" ] || mkdir "${OUTPUT_DIR}/fastas"
[ -d "${OUTPUT_DIR}/mmseqs" ] || mkdir "${OUTPUT_DIR}/mmseqs"
[ -d "${OUTPUT_DIR}/igblast" ] || mkdir "${OUTPUT_DIR}/igblast"
QUERY=$(basename "${ASSEMBLED_ANNOTATION%.*}")
TARGET=$(basename "${PROTEOMIC_DATABASE%.*}")


# Run the scripts
python3 $SCRIPTS_DIR/proteomic-database_to_cluster-rep-fasta.py \
    --proteomic_database $PROTEOMIC_DATABASE \
    --fasta ${OUTPUT_DIR}/fastas/${TARGET}.fasta
    
python3 $SCRIPTS_DIR/assembled-annotation_to_fasta.py \
    --assembled_annotation_in $ASSEMBLED_ANNOTATION \
    --fasta ${OUTPUT_DIR}/fastas/${QUERY}.fasta \
    --output_dir ${OUTPUT_DIR} \
    --assembled_annotation_out ${OUTPUT_DIR}/${QUERY}_filtered

MMSEQS_DIR="${OUTPUT_DIR}/mmseqs"
FASTA_DIR="${OUTPUT_DIR}/fastas"
mmseqs createdb $FASTA_DIR/${TARGET}.fasta $MMSEQS_DIR/${TARGET}_db
mmseqs createindex $MMSEQS_DIR/${TARGET}_db $MMSEQS_DIR/tmp
mmseqs createdb $FASTA_DIR/${QUERY}.fasta $MMSEQS_DIR/${QUERY}_db
mmseqs search $MMSEQS_DIR/${QUERY}_db $MMSEQS_DIR/${TARGET}_db $MMSEQS_DIR/${QUERY}_QUERY_${TARGET}_TARGET_db $MMSEQS_DIR/tmp --threads $THREADS
mmseqs convertalis $MMSEQS_DIR/${QUERY}_db $MMSEQS_DIR/${TARGET}_db $MMSEQS_DIR/${QUERY}_QUERY_${TARGET}_TARGET_db $OUTPUT_DIR/matches.m8 \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"

python3 $SCRIPTS_DIR/post_process.py \
    --matches $OUTPUT_DIR/matches.m8 \
    --proteomic_database $PROTEOMIC_DATABASE \
    --assembled_annotation_in ${OUTPUT_DIR}/${QUERY}_filtered \
    --peak_ratio $PEAK_RATIO \
    --output_dir $OUTPUT_DIR \
    --minimum_cluster_size 1000 \
### INCLUDE A FASTA WRITE FUNCTION IN post_process.py
    --fasta_out ${OUTPUT_DIR}/fastas/${QUERY}_filtered.fasta

${IGBLAST_PATH}/bin/igblastn \
    -germline_db_V ${IGBLAST_PATH}/database/homo_sapiens/IG_dna_clean/IGHV_clean \
    -germline_db_D ${IGBLAST_PATH}/database/homo_sapiens/IG_dna_clean/IGHD_clean \
    -germline_db_J ${IGBLAST_PATH}/database/homo_sapiens/IG_dna_clean/IGHJ_clean \
    -organism human \
    -auxiliary_data ${IGBLAST_PATH}/optional_file/human_gl.aux \
    -query ${FASTA_IN}} \
    -outfmt 19 \
    -out $OUTPUT_DIR/igblast/${FASTA_IN%.fasta}_igblast.tsv \
    -in ${OUTPUT_DIR}/fastas/${QUERY}_filtered.fasta # CHECK THIS

# CHECK HOW THIS IS RUN
python3 insertion_finder_indel_v1.py \
    --igblast_tsv $OUTPUT_DIR/igblast/${FASTA_IN%.fasta}_igblast.tsv
    #--cluster_stats $OUTPUT_DIR/cluster_stats.csv \
    #--matches $OUTPUT_DIR/matches.m8 \
    #--vh_assembled_annotation ${OUTPUT_DIR}/fastas/${QUERY}_filtered.fasta
