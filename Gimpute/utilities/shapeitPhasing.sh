#!/usr/bin/bash
#$ -cwd

CHR=$1
# parameters
THREAT=$2
EFFECTIVESIZE=$3
chrXflag=$4

# directories
ROOT_DIR=./
DATA_DIR=${ROOT_DIR}2-dataFiles/
RESULTS_DIR=${ROOT_DIR}3-phaseResults/
impRefDIR=/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/

# executable
SHAPEIT_EXEC=/data/noether/tools/imputeScript/shapeit

# reference data files
GENMAP_FILE=${impRefDIR}genetic_map_chr${CHR}_combined_b37.txt
HAPS_FILE=${impRefDIR}ALL_1000G_phase1integrated_v3_chr${CHR}_impute_macGT1.hap.gz
LEGEND_FILE=${impRefDIR}ALL_1000G_phase1integrated_v3_chr${CHR}_impute_macGT1.legend.gz
SAMPLE_FILE=${impRefDIR}ALL_1000G_phase1integrated_v3.sample

# GWAS data files in PLINK BED format
GWASDATA_BED=${DATA_DIR}gwas_data_chr${CHR}.bed
GWASDATA_BIM=${DATA_DIR}gwas_data_chr${CHR}.bim
GWASDATA_FAM=${DATA_DIR}gwas_data_chr${CHR}.fam

# main output file
OUTPUT_HAPS=${RESULTS_DIR}chr${CHR}.haps
OUTPUT_SAMPLE=${RESULTS_DIR}chr${CHR}.sample
OUTPUT_LOG=${RESULTS_DIR}chr${CHR}.log

## phase GWAS genotypes
$SHAPEIT_EXEC \
--input-bed $GWASDATA_BED $GWASDATA_BIM $GWASDATA_FAM \
--input-map $GENMAP_FILE \
$chrXflag \
--input-ref $HAPS_FILE $LEGEND_FILE $SAMPLE_FILE \
--thread $THREAT \
--effective-size $EFFECTIVESIZE \
--output-max $OUTPUT_HAPS $OUTPUT_SAMPLE \
--output-log $OUTPUT_LOG
