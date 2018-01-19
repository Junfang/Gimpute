
# Author : Junfang Chen
# Version0: 06 Jun 2016
# VersionX: 18 Jan 2018
 
 
## required functions AND libraries
source(paste0(codeDir, "4_function4imputation.R"))
library(doParallel) 

## step 1 
## Remove monomorphic SNPs from lifted/QC-ed data 

inputPrefix4aligned2impRef = "3_4_removedSnpDiffAlleles" ## will also be used in step 4 and 5;
outputPrefix = "4_1_removedMonoSnp"
outputMonoSNPfile = "4_1_snpMonoRemoved.txt" # will be used in step 4 and 5.

## copy plink files from last step; 
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".bed  ./4-imputation/") )
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".fam  ./4-imputation/") )
system( paste0("cp ./3-lifting/", inputPrefix4aligned2impRef, ".bim  ./4-imputation/") )

## remove Monomorphic SNPs
setwd('4-imputation')
removedMonoSnp(plink, inputPrefix4aligned2impRef, outputPrefix, outputMonoSNPfile)

# ## remove unwanted plink files << ## will also be used in step 4 and 5; after that you can remove.
# system(paste0('rm ', inputPrefix, '*'))


## step 2 
##########################################################################
########################################################################## imputation main pipeline

## sub-steps
# step 2.1 chrWiseSplit.R
# step 2.2 chunk4eachChr.R
# step 2.3 prePhasingByShapeit.R
# step 2.4 imputedByImpute2.R 
# step 2.5 formatConvertGtool.R 
# step 2.6 mergeImputeData.R 
# step 2.7 filterImputeData.R 
 
## One must create directories for storing temporal imputation output files 
## The name of these directories must be fixed for the sake of the subsequent steps.
tmp4imputeDIR = 'tmpImputeV4'
system( paste0("mkdir ", tmp4imputeDIR))
setwd( tmp4imputeDIR ) ## 
## sub-directories  
system( "mkdir 1-dataFiles")
system( "mkdir 2-chunkFile") 
system( "mkdir 3-phaseResults")
system( "mkdir 4-imputeResults")
system( "mkdir 5-postImpute")
system( "mkdir 6-finalResults") 
setwd("..")  


########################################
######################################## chrWiseSplit.R
## step 2.1

# copy plink files without monomorphic SNPs; prepare for the imputation.
outputPrefix = "4_1_removedMonoSnp"
outputPrefixTmp  = 'gwas_data_chr'  

system( paste0("scp ", outputPrefix, ".bed  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".bed"))
system( paste0("scp ", outputPrefix, ".fam  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".fam"))
system( paste0("scp ", outputPrefix, ".bim  ./", tmp4imputeDIR, "/1-dataFiles/", outputPrefixTmp, ".bim"))
  

setwd(paste0('./', tmp4imputeDIR, '/1-dataFiles/'))
inputPrefix  = outputPrefixTmp 
outputPrefix  = outputPrefixTmp ## Same as the input files; but the chromosome number will be appended
chrX_PAR1suffix = 'X_PAR1'
chrX_PAR2suffix = 'X_PAR2'
PAR = chrWiseSplit(plink, inputPrefix, chrX_PAR1suffix, chrX_PAR2suffix)
 
PAR 
 
if (PAR[[1]]) {par1 = 'X_PAR1'} else {par1 = NULL}
if (PAR[[2]]) {par2 = 'X_PAR2'} else {par2 = NULL}

  

########################################
######################################## chunk4eachChr.R
## step 2.2

inputPrefix = 'gwas_data_chr'
outputPrefix = 'chunks_chr' 
chrs = c(1:23, par1, par2)   
windowSize = 3000000 
chunk4eachChr(inputPrefix, outputPrefix, chrs, windowSize) 

setwd('..') 
system( paste0("mv ./1-dataFiles/", outputPrefix, "*.txt  ./2-chunkFile/"))



## step 2.3  

# define directories
dataDIR = "./1-dataFiles/"
impRefDIR = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
phaseDIR = "./3-phaseResults/"
 
prefix4plinkEachChr = 'gwas_data_chr'  
# chrs = c(1, 11, 23, par1, par2)   
nThread = 40
effectiveSize = 20000 
nCore = 1 
prePhasingByShapeit(shapeit, chrs, dataDIR, prefix4plinkEachChr, impRefDIR, phaseDIR, nThread, effectiveSize, nCore)

 

# step 2.4   
# define directories
prefixChunk = "./2-chunkFile/chunks_chr"  ## 
phaseDIR = "./3-phaseResults/"
impRefDIR = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/"
imputedDIR = "./4-imputeResults/"

prefix4plinkEachChr = 'gwas_data_chr'
# chrs = c(1,23, par1, par2)     
nCore = 30 ## try to tune it for your own data size
effectiveSize = 20000  

imputedByImpute2(impute2, chrs, prefixChunk, phaseDIR, impRefDIR, imputedDIR, prefix4plinkEachChr, nCore, effectiveSize)


  
 
# step 2.5   
####################################################### extract only SNPs (without INDELs)
setwd("./4-imputeResults")  
## extract only SNPs starting with 'rs';  .
ls = system("ls gwas*.impute2", intern=T)
variantPrefix = 'rs' 

biglists = as.list(ls)
mclapply(biglists, function(i){
	arg1 = paste0(i, "noINDEL.impute2")
	arg2 = paste0("grep '", variantPrefix, "' ", i, " | awk '{if(length($4)==1 && length($5)==1) print}' > ", arg1)
	# print(arg2)
	system(arg2)
}, mc.cores= 38) 
setwd("..")  
####################################################### <<<


# chrs = c(1,23, par1, par2)
   
prefixChunk = "./2-chunkFile/chunks_chr"  ## 
phaseDIR = "./3-phaseResults/" 
imputedDIR = "./4-imputeResults/"
prefix4plinkEachChr = 'gwas_data_chr'
suffix4imputed = ".impute2noINDEL.impute2"
# suffix4imputed = ".impute2"
postImputeDIR = './5-postImpute/' 

nCore = 30   
formatConvertGtool(gtool, chrs, prefixChunk, phaseDIR, imputedDIR, prefix4plinkEachChr, suffix4imputed, postImputeDIR, nCore)

 


# step 2.6  
####################################################### Modify missing genotype format.
setwd("./5-postImpute/")  
# replace 'N' in the .ped files into 0 > missing values.
chrslist = as.list(chrs)
prefix4plinkEachChr = 'gwas_data_chr' ## just for parallel computing
fn = mclapply(chrslist, function(i){
	system(paste0("sed -i 's/N/0/g' ", prefix4plinkEachChr, i, ".*ped "))
}, mc.cores=nCore)
## Note: check if also any "N" in *.fam files. If so, change back after merging. 
####################################################### <<< 

prefix4plinkEachChr = 'gwas_data_chr'
prefix4imputedPlink = 'gwasImputed'
par1 = 'X_PAR1'
par2 = 'X_PAR2'
# chrs = c(1,23, par1, par2)   
mergeImputeData(plink, chrs, prefix4plinkEachChr, prefix4imputedPlink, nCore)
 
 

#######################################################
setwd("..")
system("mv ./5-postImpute/gwasImputed* ./6-finalResults/") 
system("mv ./4-imputeResults/*.impute2_info ./6-finalResults/") 

# step 2.7    
setwd("./6-finalResults")
suffix4impute2info = ".impute2_info"
impute2infoFile = "impute2infoUpdated.txt"
infoScore = 0.6
badImputeSNPfile = "badImputeSNPs.txt"
prefix4imputedPlink = "gwasImputed"
prefix4imputedFilterPlink = "gwasImputedFiltered" 
filterImputeData(plink, suffix4impute2info, impute2infoFile, infoScore, badImputeSNPfile, prefix4imputedPlink, prefix4imputedFilterPlink)

setwd("..")
setwd("..")

########################################################################## After imputation

## step 2 
## Final imputed results >> 
 
prefix4imputedPlink = "gwasImputed"
## output file name change to: 
imputedDatasetfn = "4_2_imputedDataset"
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".bed ", imputedDatasetfn, ".bed") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".bim ", imputedDatasetfn, ".bim") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedPlink, ".fam ", imputedDatasetfn, ".fam") )


## step 3
## Filtered imputed data set; Remove imputed SNPs with (info < 0.6), only retain 'Good' SNPs.
prefix4imputedFilterPlink = "gwasImputedFiltered"
filteredImputedDatasetfn = "4_3_removedSnpInfoPostImp" 
snpWithBadInfoFile = "4_3_snpRemovedInfoPostImp.txt"
snpImputedInfoScoreFile = "4_3_snpImputedInfoScore.txt"
 
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".bed ", filteredImputedDatasetfn, ".bed") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".bim ", filteredImputedDatasetfn, ".bim") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", prefix4imputedFilterPlink, ".fam ", filteredImputedDatasetfn, ".fam") )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", impute2infoFile, " ", snpImputedInfoScoreFile) )
system( paste0("scp ./", tmp4imputeDIR, "/6-finalResults/", badImputeSNPfile, " ", snpWithBadInfoFile) )


## step 4 
## Remove previous identified monomorphic SNPs in the imputed dataset.
filteredImputedDatasetfn = "4_3_removedSnpInfoPostImp" 
removedMonoSnpAfter = "4_4_removedMonoSnpAfter"
# perl4snpsWithSamePos = "/home/junfang.chen/groupSVN/trunk/datasets/bin/dataProcess/4-imputation/snpsWithSamePos.pl"

## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", removedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", removedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", removedMonoSnpAfter, ".fam") )  
} else { 
	## extract PLINK files contain only monomorphic SNPs from the original aligned (lifted and QC-ed) data set.
	system( paste0(plink, " --bfile ", inputPrefix4aligned2impRef, " --extract ", outputMonoSNPfile, " --make-bed --out ", inputPrefix4aligned2impRef, "Tmp") ) 
	system( paste0(perl4snpsWithSamePos, " ", filteredImputedDatasetfn, ".bim ", inputPrefix4aligned2impRef, "Tmp", ".bim > tmp.txt" ))
	system( paste0(plink, " --bfile ", filteredImputedDatasetfn, " --exclude tmp.txt --make-bed --out ", removedMonoSnpAfter) )
} 


## step 5
## Add previous identified monomorphic SNPs in the imputed dataset.
addedMonoSnpAfter = "4_5_addedMonoSnpAfter" 

 ## if no monomorphic SNPs:
if ( file.size(paste0(outputMonoSNPfile))==0 ){ 
		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", addedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", addedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", addedMonoSnpAfter, ".fam") )  
} else { 
  
	## merge both datasets
	system( paste0(plink, " --bfile ", removedMonoSnpAfter, " --bmerge ", 
		inputPrefix4aligned2impRef, ".bed ", inputPrefix4aligned2impRef, ".bim ", inputPrefix4aligned2impRef, ".fam --make-bed --out ", addedMonoSnpAfter) )
	## remove tmp files
	# system( paste0("rm tmp.txt") )
	# system( paste0("rm ", inputPrefix4aligned2impRef, "*") )
}  


## step 6
## Remove SNPs which have a non missing value for less then 20 instances. 

inputPrefix = addedMonoSnpAfter  
missCutoff = 20
outputPrefix = "4_6_removedSnpMissPostImp"
snpWithManyMissSNPfile = "4_6_snpRemovedMissPostImp.txt"
removedSnpMissPostImp(plink, inputPrefix, missCutoff, snpWithManyMissSNPfile, outputPrefix)

   
setwd("..")
  