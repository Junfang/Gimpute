
  
# File   : 3_lifting.R
# Author : Junfang Chen
# Version0: 06 Jun 2016
# VersionX: 15 Jan 2018
 



## required functions 
source(paste0(codeDir, "4_function4imputation.R"))



## create directories for storing imputation output files
setwd("./4-imputation")
system( "mkdir tmp4impute")
setwd( "./tmp4impute") ## 
## sub-directories
system( "mkdir 1-checkAlign")
system( "mkdir 2-dataFiles")
system( "mkdir 3-phaseResults")
system( "mkdir 4-chunkFile")
system( "mkdir 5-imputeResults")
system( "mkdir 6-postImpute")
system( "mkdir 7-finalResults")
setwd("..")  ## go to ./4-imputation 
setwd("..")  ## go to main directory 

 

## step 1 
## Remove monomorphic SNPs from lifted/QC-ed data 
inputPrefix = "3_4_removedSnpDiffAlleles"
outputPrefix = "4_1_removedMonoSnp"
outputPrefix.snp = "4_1_snpMonoRemoved"
outputPrefixTmp = "gwas_data_chr"
removedMonoSnp(plink, inputPrefix, outputPrefix, outputPrefix.snp, outputPrefixTmp)
 
    
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


## step 2.1
inputPrefix  = 'gwas_data_chr'  
X_PAR1sufix = 'X_PAR1'
X_PAR2sufix = 'X_PAR2'
par = chrWiseSplit(plink, inputPrefix, X_PAR1sufix, X_PAR2sufix)
par 
 
if (par[[1]]) {par1 = 'X_PAR1'} else {par1 = NULL}
if (par[[2]]) {par2 = 'X_PAR2'} else {par2 = NULL}


## step 2.2
inputPrefix = 'gwas_data_chr'
outputPrefix = 'chunks_chr' 
chrs = c(21:23, par1, par2)   
windowSize = 3000000 
chunk4eachChr(inputPrefix, outputPrefix, chrs, windowSize) 
 

## step 2.3
## preshasing scripts
prePhasingScript = "/home/junfang.chen/Gimpute/utilities/shapeitPhasing.sh"
##      
chrs = c(21:23, par1, par2)   
nThreat = 40
effectiveSize = 20000 
chrXflag = "--chrX"
prePhasingByShapeit(prePhasingScript, chrs, nThreat, effectiveSize, chrXflag) 



# step 2.4   
## imputation script
impute2script = "/home/junfang.chen/Gimpute/utilities/impute2.sh"
impute2script4chrX = "/home/junfang.chen/Gimpute/utilities/impute2chrX.sh"
##      
chrs = c(21:23, par1, par2)    
chrs = c(21:22, par1, par2)    

nCore = 30
effectiveSize = 20000 
XPAR = "-Xpar"
imputedByImpute2(impute2script, chrs, effectiveSize, nCore, XPAR) 
 


# step 2.5   
## conversion scripts 
gtoolScript = "/home/junfang.chen/Gimpute/utilities/gtoolOutput.sh"
##      
variantPrefix = "rs"
chrs = c(21:23, par1, par2)   
nCore = 30 
formatConvertGtool(variantPrefix, chrs, gtoolScript, nCore)

 
# step 2.6  
chrs = c(1:23, par1, par2)   
mergeImputeData(plink, chrs)


# step 2.7    
filterImputeData(plink, infoScore=0.6)


## step 3
##########################################################################
########################################################################## Data re-organization  

## define the file name generated after imputation. 
snpMonoRemoved = "4_1_removedMonoSnpBefore" 
snpMonoRemoved.snp = "4_1_snpMonoRemoved" 
 

imputedDatasetfn = "4_2_imputedDataset"
filteredImputedDatasetfn = "4_3_removedSnpInfoPostImp"
snpWithBadInfo = "4_3_snpRemovedInfoPostImp"
 

# 4_4_removedMonoSnpAfter.%
removedMonoSnpAfter = "4_4_removedMonoSnpAfter"
# 4_5_addedMonoSnpAfter.% 
addedMonoSnpAfter = "4_5_addedMonoSnpAfter"
 
removedSnpMissPostImpfn = "4_6_removedSnpMissPostImp"
snpWithManyMissSNPsfn = "4_6_snpRemovedMissPostImp"
## /4-imputation
## 4.2 Final imputed results (fam file has to be updated)
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputed.bed ", imputedDatasetfn, "_tmp.bed") )
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputed.bim ", imputedDatasetfn, "_tmp.bim") )
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputed.fam ", imputedDatasetfn, "_tmp.fam") )

# In some cases, one has to update the IDs in the .fam files after imputation.
## update sample IDs in the imputed fam file by using original fam file information
## changes ID codes for individuals specified in recoded.txt, 
## which should be in the format of four columnds per row: old FID, old IID, new FID, new IID, e.g. 
famNew = read.table(paste0(snpMonoRemoved, '.fam'), stringsAsFactors=F)
head(famNew) 
famOld = read.table(paste0(imputedDatasetfn, '_tmp.fam'), stringsAsFactors=F)
head(famOld) 
fam2update = cbind(famOld[,1:2], famNew[,1:2])
print(sum(famOld[,1]==famNew[,2]))

write.table(fam2update, file="famIDupdate.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
system(paste0(plinkB,  imputedDatasetfn, "_tmp", " --update-ids famIDupdate.txt --make-bed --out ", imputedDatasetfn  )) 
system(paste0('rm ', imputedDatasetfn, "_tmp*"))






## 4.3 Filtered imputed results (fam file has to be updated)
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputedFiltered.bed ", filteredImputedDatasetfn, "_tmp.bed") )
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputedFiltered.bim ", filteredImputedDatasetfn, "_tmp.bim") )
system( paste0("scp ./tmp4impute/7-finalResults/gwasImputedFiltered.fam ", filteredImputedDatasetfn, "_tmp.fam") )

## update sample IDs in the imputed fam file by using original fam file information
system(paste0(plinkB,  filteredImputedDatasetfn, "_tmp", " --update-ids famIDupdate.txt --make-bed --out ", filteredImputedDatasetfn  )) 

system(paste0('rm ', filteredImputedDatasetfn, "_tmp*"))
system(paste0("rm famIDupdate.txt"))
 





# 4. Remove the monomorphic SNPs in 4_1_snpMonoRemoved.txt.
plinkB = "/data/noether/tools/plink/plink --bfile " 

if ( file.size(paste0(snpMonoRemoved.snp, ".txt"))==0 ){ 
		
		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", removedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", removedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", removedMonoSnpAfter, ".fam") ) 

		system(paste0("scp ", filteredImputedDatasetfn, ".bed ", addedMonoSnpAfter, ".bed") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".bim ", addedMonoSnpAfter, ".bim") ) 
		system(paste0("scp ", filteredImputedDatasetfn, ".fam ", addedMonoSnpAfter, ".fam") ) 

} else { 


	setwd("..")
	inputLifted = "3_4_removedSnpDiffAlleles"

	system( paste0("scp ./3-lifting/", inputLifted, ".bed  ./4-imputation/") )
	system( paste0("scp ./3-lifting/", inputLifted, ".bim  ./4-imputation/") )
	system( paste0("scp ./3-lifting/", inputLifted, ".fam  ./4-imputation/") )
	setwd("./4-imputation")
 
	## generate/extract mono SNPs from the raw lifted data
	system( paste0(plinkB, inputLifted, " --extract ", snpMonoRemoved.snp, ".txt --make-bed --out ", inputLifted, "Tmp") ) 

	 
	dir4pl = "/home/junfang.chen/groupSVN/trunk/datasets/bin/dataProcess/4-imputation/"
	system( paste0(dir4pl, "snpsWithSamePos.pl ", filteredImputedDatasetfn, ".bim ", inputLifted, "Tmp", ".bim > tmptxt" ))
	system( paste0(plinkB, filteredImputedDatasetfn, " --exclude tmp.txt --make-bed --out ", removedMonoSnpAfter) )



	# 5. Add the monomorphic SNPs in 4_1_snpMonoRemoved.txt with their values from the lifted instance data.
	## copy the momo genotypes from lifted data 
	system( paste0(plinkB, removedMonoSnpAfter, " --bmerge ", inputLifted, ".bed ", inputLifted, ".bim ", inputLifted, ".fam --make-bed --out ", addedMonoSnpAfter) )
	## remove tmp files
	system( paste0("rm tmp", ".txt") )
	system( paste0("rm ", inputLifted, "*") )

} 

# 6. Remove SNPs which have a non missing value for less then 20 instances.
# Primary output files: 4_4_removedSnpMissPostImp.#
# Output file with the removed SNP names: 4_4_snpRemovedMissPostImp.txt

## get the missing info 
system(paste0(plinkB, addedMonoSnpAfter, " --missing --out ", addedMonoSnpAfter)) 

missSNPinfo = read.table(paste0(addedMonoSnpAfter, ".lmiss"), stringsAsFactors=F, h=T)
# dim(missSNPinfo)
# head(missSNPinfo)
missSNPinfo[,6] <- missSNPinfo[,"N_GENO"] - missSNPinfo[,"N_MISS"] 
snpWithManyMissSNPs <- missSNPinfo[which(missSNPinfo[,6] < 20), "SNP"]
# str(snpWithManyMissSNPs) 
write.table(snpWithManyMissSNPs, file=paste0(snpWithManyMissSNPsfn, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")

system( paste0(plinkB, addedMonoSnpAfter, " --exclude ", snpWithManyMissSNPsfn, ".txt --make-bed --out ", removedSnpMissPostImpfn) )
system( "rm *.imiss *.lmiss *.log") 

# 
 
## make new directories 
setwd("..")

# system('mkdir 5-reductAndExpand')
# system('mkdir 6-finalResults')

 