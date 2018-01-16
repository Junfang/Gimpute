  
  
# File   : 2_qualityControl.R
# Author : Junfang Chen
# Version0: 08 Jun 2016
# VersionX: 08 Jan 2018

 
  
## required functions
source(paste0(codeDir, "2_function4QC.R"))



## step 0
## copy the last output plink files from 1-conversion
inputPrefix4QC = "1_11_removedYMtSnp"
system( paste0("scp ./1-conversion/", inputPrefix4QC, ".*", " ./2-QC/") )
setwd("./2-QC/")

## step 1
inputPrefix = inputPrefix4QC
cutoff = 0.005
outputPrefix = "2_01_removedSnpHetX"
outputPrefix.snp = "2_01_snpHetXInstNumber" 
removedSnpHetX(plink, inputPrefix, cutoff, outputPrefix, outputPrefix.snp)


## step 2  2_02_removedHetXInst
inputPrefix = "2_01_removedSnpHetX"
cutoff = 15
outputPrefix = "2_02_removedInstHetX"
outputPrefix.snp = "2_02_instHetXSnpNumber" 
removedInstHetX(plink, inputPrefix, cutoff, outputPrefix, outputPrefix.snp)
 

## step 3 
# 3. Set all heterozygous alleles of SNPs of the chromosome 23 for males
inputPrefix = "2_02_removedInstHetX" 
outputPrefix = "2_03_setHeteroHaploMissing" 
setHeteroHaploMissing(plink, inputPrefix, outputPrefix)
 

## step 4   SNP missingness < 0.05 (before sample removal);  
inputPrefix = "2_03_setHeteroHaploMissing" 
snpMissing = 0.05 #
outputPrefix = "2_04_removedSnpMissPre" 
removedSnpMissPre(plink, snpMissing, inputPrefix, outputPrefix)

## step 5 
# subject missingness < 0.02; 
inputPrefix = "2_04_removedSnpMissPre" 
sampleMissing = 0.02
outputPrefix = "2_05_removedInstMiss" 
removedInstMiss(plink, sampleMissing, inputPrefix, outputPrefix)
 
## step 6 
inputPrefix = "2_05_removedInstMiss" 
Fhet = 0.2 #
outputPrefix = "2_06_removedInstFhet" 
removedInstFhet(plink, Fhet, inputPrefix, outputPrefix)
  	

## step 7 
inputPrefix = "2_06_removedInstFhet" 
kinshipValue = 0.11 #
outputPrefix = "2_07_removedInstRelated" 
outputPrefix.ID = "2_07_instRemovedKinships" 
removedInstRelated(plink, kinshipValue, inputPrefix, outputPrefix, outputPrefix.ID)


## step 8 
inputPrefix = "2_07_removedInstRelated"  
outputPrefix = "2_08_removedParentIdsMiss" 
removedParentIdsMiss(inputPrefix, outputPrefix)


## step 9 
snpMissing2 = 0.02
inputPrefix = "2_08_removedParentIdsMiss"  
outputPrefix = "2_09_removedSnpMissPost" 
removedSnpMissPost(plink, snpMissing2, inputPrefix, outputPrefix)


## step 10 
 ##  Remove SNPs with difference >= 0.02 of SNP missingness between cases and controls.
inputPrefix = "2_09_removedSnpMissPost"  
outputPrefix = "2_10_removedSnpMissDiff" 
removedSnpMissDiff(plink, snpMissing2, inputPrefix, outputPrefix)

## step 11 
femaleChrXmiss = TRUE
inputPrefix = "2_10_removedSnpMissDiff"  
outputPrefix = "2_11_removedSnpMissAllo" 
removedSnpFemaleChrXmiss(plink, femaleChrXmiss, inputPrefix, outputPrefix)
  

## step 12
hweAUTOcontrol = TRUE
pval = 0.000001
inputPrefix = "2_11_removedSnpMissAllo"  
outputPrefix = "2_12_removedSnpHweAutoCt" 
outputPrefix.pval = "2_12_snpHweValuesAutoCt" 
outputPrefix.snp = "2_12_snpRemovedHweAutoCt" 
removedSnpHWEautoControl(plink, hweAUTOcontrol, pval, inputPrefix, outputPrefix)



## step 13
femaleChrXhweControl = TRUE #  
# femaleChrXhweAll = TRUE #   
pval = 0.000001
inputPrefix = "2_12_removedSnpHweAutoCt"  
outputPrefix = "2_13_removedSnpHweAlloCt" 
outputPrefix.pval = "2_13_snpHweValuesAlloCt" 
outputPrefix.snp = "2_13_snpRemovedHweAlloCt" 
removedSnpFemaleChrXhweControl(plink, femaleChrXhweControl, pval, inputPrefix, outputPrefix)




 

## in order to get the ethnic group info
setwd('..')
metaDataFile = "1_01_metaData.txt"
system( paste0("scp ./1-conversion/", metaDataFile, " ./2-QC/") )
setwd("./2-QC/")


## step 14 

## 13  Remove chrX SNPs with HWE (p < 10.6) in female controls.
inputPrefix = "2_13_removedSnpHweAlloCt"
outputPCs = "2_14_eigenvecQCedBeforeRm"
plotPCA4plink(gcta, inputPrefix, outputPCs)

## remove outliers
inputPCs = "2_14_eigenvecQCedBeforeRm"
cutoff =  NULL ##  
cutoffSign = 'smaller'
outputIDs = "2_14_eigenvecQCedRemoved"
inputPrefix = "2_13_removedSnpHweAlloCt"
outputPCplot = "2_14_eigenvecQCedAfterRm"
outputPrefix = "2_14_removedOutliersQCed" 
removeOutlierByPCs(inputPCs, cutoff, cutoffSign, outputIDs, plink, gcta, inputPlink, outputPCplot, outputPrefix)
 
 
######################## 
## remove unwanted files
## remove meta file and metaDataFile+AA.txt
system( paste0("rm  ", metaDataFile) ) 
system( paste0("rm  *.log *.hh") )
## change dir 
setwd("..") 


  
