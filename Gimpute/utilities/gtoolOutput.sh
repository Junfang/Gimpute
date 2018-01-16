#!/usr/bin/bash
#$ -cwd
 
CHR=$1
CHUNK_START=`printf "%.0f" $2`
CHUNK_END=`printf "%.0f" $3`
 
# directories

ROOT_DIR=./

PhaseRESULTS_DIR=${ROOT_DIR}3-phaseResults/
imputeRESULTS_DIR=${ROOT_DIR}5-imputeResults/
OUTPUT_DIR=${ROOT_DIR}6-postImpute/

# executable (must be version 2.2.0 or later)
GTOOL_EXEC=/data/noether/tools/imputeScript/gtool
 
 
## MODIFY THE FOLLOWING THREE LINES TO ACCOMODATE OTHER PANELS
# INPUT data files
SAM_FILE=${PhaseRESULTS_DIR}chr${CHR}.sample
GEN_FILE=${imputeRESULTS_DIR}gwas_data_chr${CHR}.pos${CHUNK_START}-${CHUNK_END}.impute2noINDEL.impute2


# METHOD-B: haplotypes from SHAPEIT phasing run
# GWAS_HAPS_FILE=${RESULTS_DIR}chr${CHR}.haps
 
# main output file
PED_FILE=${OUTPUT_DIR}gwas_data_chr${CHR}.pos${CHUNK_START}-${CHUNK_END}.ped
MAP_FILE=${OUTPUT_DIR}gwas_data_chr${CHR}.pos${CHUNK_START}-${CHUNK_END}.map

# OUTPUT_FILE=${OUTPUT_DIR}gwas_data_chr${CHR}.pos${CHUNK_START}-${CHUNK_END}.gtool
 
## impute genotypes from GWAS haplotypes
$GTOOL_EXEC -G \
   --g $GEN_FILE \
   --s $SAM_FILE \
   --phenotype plink_pheno \
   --chr ${CHR} \
   --ped $PED_FILE \
   --map $MAP_FILE \
   --snp
   # --log $OUTPUT_FILE \