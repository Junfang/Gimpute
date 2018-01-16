
# File   : 0_dataImputePipeline.R
# Author : Junfang Chen
# Version0: 28 Jun 2016
# VersionX: 11 Jan 2018

## This script is used to run the whole genetic data imputation process including 
## data cleaning, snp information update (lifting), QC, outlier detection,  
## alignment check, Imputation and post imputation.
 



############################################################
### code chunk number 0: create all necessary directories   
############################################################
## 
system('mkdir 0-rawData')
system('mkdir 1-conversion')
system('mkdir 2-QC')
system('mkdir 3-lifting')
system('mkdir 4-imputation')
system('mkdir 5-reductAndExpand')
system('mkdir 6-finalResults')

system('cd 0-rawData')
system('mkdir plinkFiles') ## Original plink files 
system('mkdir sampleInfo') ## Meta data information
system('cd ..')


## required global libraries/tool paths/parameters


## plink path
plink = "/data/noether/tools/plink/plink "  
gcta = "/data/noether/tools/gcta/gcta64 " 

## the script directory 
codeDir = "/home/junfang.chen/Gimpute/R/"

 

############################################################
### code chunk number 1: SNP information update  
############################################################
## processed data will be generated in this directory: ./1-conversion/
source("1_snpInfoUpdate.R")


############################################################
### code chunk number 2: Quality Control
############################################################
## processed data will be generated in this directory: ./2-QC/
source("2_qualityControl.R")



############################################################
### code chunk number 3: Quality control   
############################################################
source("3_lifting.R")



############################################################
### code chunk number 4: Imputation
############################################################
source("4_imputation.R")


############################################################
### code chunk number 5: Data subset and expansion 
############################################################
source("5_reductAndExpand.R")


############################################################
### code chunk number 6: Final result
############################################################
source("6_outputFinal.R")
 
 