 
# File   : snpInfoUpdate.R
# Author : Junfang Chen
# Version0: 08 Jun 2016
# VersionX: 11 Jan 2018

# This script is used for updating snp information from one genome build to another.
# for the Affymetric and Illumina genotype data

 
## required tool: plink 


## required functions
source(paste0(codeDir, "1_function4infoUpdate.R"))

## pre-defined parameters
# rawPlinkFiles: the prefix of input plink file name
rawPlinkFiles = 'snp270'
dupSampleIDs = NULL
excludedProbeIds = NULL  
chipAnnoFile= "/home/junfang.chen/Gimpute/config/Illumina/Human1M-Duov3_B-b37.Illmn.strand"
 

########################### go to main directory
 
## step 1
## copy plink files and meta information file
system( paste0("scp ./0-rawData/plinkFiles/", rawPlinkFiles, ".*  ./1-conversion/") ) 
system( paste0("scp ./0-rawData/sampleInfo/1_01_metaData.txt ./1-conversion/") ) 
setwd("./1-conversion/") 

##############################################################
############################################################## snpInfoUpdate  

## step 2  
dupSampleIDs.svn = dupSampleIDs
inputPrefix = rawPlinkFiles
outputPrefix = "1_02_removedExclInst" 
removeDupID(plink, dupSampleIDs.svn, inputPrefix, outputPrefix)
 

# step 3 replace group IDs
 
metaDataFile = "1_01_metaData.txt"
inputPrefix = "1_02_removedExclInst"
outputPrefix = "1_03_replacedGroupAndSex"
replaceGroupId(plink, inputPrefix, metaDataFile, outputPrefix)

# step 4 remove instances without group IDs
metaDataFile = "1_01_metaData.txt"
inputPrefix = "1_03_replacedGroupAndSex"
outputPrefix = "1_04_removedNoGroupId"
removeNoGroupId(plink, inputPrefix, outputPrefix)

## step5 remove instances with improper ancestry 
metaDataFile = "1_01_metaData.txt"
ancestrySymbol = 'EA'
inputPrefix = "1_04_removedNoGroupId"
outputPrefix = "1_05_removedWrongAnceInst"
removedWrongAnceInst(plink, inputPrefix, metaDataFile, ancestrySymbol, outputPrefix)

## step 6 
inputPrefix = "1_05_removedWrongAnceInst"
excludedProbeIds = excludedProbeIds  ## excludedProbeIds is defined in svn
outputPrefix = "1_06_removedExclProbe" ## 
removedExclProbe(plink, inputPrefix, excludedProbeIds, outputPrefix) 


## step 7 
inputPrefix = "1_06_removedExclProbe"
chipAnnoFile.svn = chipAnnoFile ## defined in SVN
chipType = 'illumina'
outputPrefix = "1_07_removedUnmapProbes"   
outputPrefix.snp = "1_07_probesUnmapped2ChipRef"

removedUnmapProbes(plink, inputPrefix, chipAnnoFile.svn, outputPrefix, outputPrefix.snp)

## step 8 
inputPrefix = "1_07_removedUnmapProbes" 
chipAnnoFile.svn = chipAnnoFile ## defined in SVN
chipType = 'illumina'
outputPrefix = "1_08_removedDoubleProbes"   
outputPrefix.snp = "1_08_probesDouble"
removedDoubleProbes(plink, inputPrefix, chipAnnoFile.svn, chipType, outputPrefix, outputPrefix.snp)
 
## step 9
inputPrefix = "1_08_removedDoubleProbes" 
chipAnnoFile.svn = chipAnnoFile ## defined in SVN
chipType = 'illumina'
inputPrefix.rs = "1_09.updatedSnp2rs" ## tobeRemoved
inputPrefix.chr = "1_09.updatedSnpchr" ## tobeRemoved
inputPrefix.pos = "1_09.updatedSnppos" ## tobeRemoved
inputPrefix.strand = "1_09.updatedSnpstrand" ## tobeRemoved
outputPrefix = "1_09_updatedSnpInfo"   
updatedSnpInfo(plink, inputPrefix,  chipAnnoFile.svn, chipType, inputPrefix.rs, inputPrefix.chr, inputPrefix.pos, inputPrefix.strand, outputPrefix)

## step 10 
inputPrefix = "1_09_updatedSnpInfo"
outputPrefix = "1_10_changedXyChr"
changedXyChr(plink, inputPrefix, outputPrefix)

## step 11 
inputPrefix = "1_10_changedXyChr"
outputPrefix = "1_11_removedYMtSnp"
removedYMtSnp(plink, inputPrefix, outputPrefix)


## remove unwanted files
system( paste0("rm  *.log") )
system( paste0("rm  *.hh") )
## change dir to the main directory
setwd("..")   










