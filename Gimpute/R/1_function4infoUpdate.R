 
# File   : function4infoUpdate.R
# Author : Junfang Chen
# Version0: 08 Jun 2016
# VersionX: 08 Jan 2018

## This script only contains r functions that are used to update/convert/lift genotypes 
## from one genome build to  another including basic QC: data cleaning. 
## Note the difference between chip types.
 

#########################################################################################  
######################################################################################### 
## snp information conversion/updating/re-mapping

##########################################   
########################################## replace chr name
## If the plink files contain names of the chromosomes X, Y, 
## XY (pseudo-autosomal region of X) and MT (mitochondrial) 
## replace these by the following numbers: 
## X by 23, Y by 24, XY by 25 and MT by 26.
## plink: --update-chr [filename] {chr col. number} {variant ID col.}
renamedSnpChr <- function(plink, inputPrefix, outputPrefix){

	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=F)
	whX = which(bim[,1]=="X")
	whY = which(bim[,1]=="Y")
	whXY = which(bim[,1]=="XY")
	whMT = which(bim[,1]=="MT")

 	tmpBIM = bim
 	tmpBIM[whX,1] <- 23
 	tmpBIM[whY,1] <- 24
 	tmpBIM[whXY,1] <- 25
 	tmpBIM[whMT,1] <- 26
 	
 	updateSNPchr = tmpBIM[c(whX, whY, whXY, whMT), c(2,1)] ## change to --> 1st SNP; 2nd chr names
  	write.table(updateSNPchr, file=paste0(outputPrefix, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
 	system( paste0(plink, " --bfile ", inputPrefix,  " --update-chr ", paste0(outputPrefix, ".txt"), " 2 1 --make-bed --out ", outputPrefix) )  
 	# system( paste0("rm ", outputPrefix, ".txt") )  
}

 
 

##########################################   
##########################################  removeDupID.R
#### clean/prepare the raw plink file  
# in the FAM file
# colnames(fam) = c("famid", "id", "paternalID", "maternalID", "sex", "phenotype") 
# check the followings
# 1. any .dup IDs (both famid and id), if so remove IDs with .dup
# 2. any duplicated IDs (only for id), if so remove them
# 3. correct sex format?  (1=male; 2=female; other=unknown) 
# 4. correct phenotype ? Assuming a disease phenotype (1=unaff, 2=aff, 0=miss) 

## remove .dup individuals 
 
removeDupID <- function(plink, dupSampleIDs.svn, inputPrefix, outputPrefix){

	if (!is.null(dupSampleIDs.svn)){

		famv0 = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE)  
		dupIDs = read.table(file=dupSampleIDs.svn, stringsAsFactors=FALSE)  
		## write into plink format .txt for removal
		dupIDs4plink = famv0[is.element(famv0[,2], dupIDs[,1]), 1:2] 
		dupIDsfn = paste0(outputPrefix, ".txt")
		write.table(dupIDs4plink, file=dupIDsfn, quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 

		system( paste0(plink, " --bfile ", inputPrefix, " --remove ", dupIDsfn, " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", dupSampleIDs.svn) )
		system( paste0("rm ", dupIDsfn) )

		## remove the raw genotype plink files
		system( paste0("rm ", inputPrefix, ".*") )

	} else { 
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

		## remove the raw genotype plink files
		system( paste0("rm ", inputPrefix, ".*") )
	  
	    } 
}

 
##########################################   
##########################################  replaceGroupId.R


## 1. find the shared Ids between plink files and metadata file 
##    a. update the GroupId to the right format
## 	  b. Group should be 1 and 2. (1=unaff, 2=aff, 0=miss) instead of case 1 control 0
# metaData: data.frame with at least 2 columns, 1st column is sample IDs, 2nd is the group, case 1 control 0.
  

## 1. check the phenotypes, if there is any sample without phenotype, then remove it.
## for some datasets, there is no such case
# metaData: data.frame with at least 2 columns, 1st column is sample IDs, 2nd is the group, case 1 control 0.
 

replaceGroupId <- function(plink, inputPrefix, metaDataFile, outputPrefix){
 
	fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
	metaData = read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
	## to make sure only compare with the same set of plIndID --> IID
	interIDs = intersect(fam[,2], metaData[,"IID"]) 	# the shared IDs 
	metadataSub = metaData[is.element(metaData[,"IID"], interIDs),]

	## one must keep the order of .fam unchanged!!
	metadataSubsort = metadataSub[match(fam[,2], metadataSub[,"IID"]),]
	missIDsIndex = which(is.na(metadataSubsort[,1])==TRUE) ## IID is in the 1st col; 
	fam[,6] = metadataSubsort[,"group"] + 1
	fam[missIDsIndex, 6] = -9  ## keep the IDs of no missing group info as -9
 	newPheno = fam[,c(1,2,6)]  ##  a pheno file that contains 3 columns (one row per individual)
 	write.table(newPheno, file="pheno.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")  

	## Alternate phenotype files
	system( paste0(plink, " --bfile ", inputPrefix, " --pheno pheno.txt --make-bed --out ", outputPrefix) )
	system( "rm pheno.txt" )

}

 
##########################################   
########################################## removeNoGroupId.R

removeNoGroupId <- function(plink, inputPrefix, outputPrefix){
 
	fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
	## 1. check if any sample without phenotypes
	noGroupIds = fam[which(fam[,6] == -9), 1:2]
	print(dim(noGroupIds))
	## if any then remove and afterwards add phenotypes
	noGroupIdsfn = paste0(outputPrefix, ".txt")
	write.table(noGroupIds, file=noGroupIdsfn, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")   
	system( paste0(plink, " --bfile ", inputPrefix, " --remove ", noGroupIdsfn, " --make-bed --out ", outputPrefix) )  

	system( paste0("rm ",  noGroupIdsfn) )
}
 
  
 
## step 5
## wrong or improper ancestry instances define in the meta data >> 
## ancestrySymbol == 'EA' In this case, only keep EA 
## If you want to impute your data set using Multi-Population Reference Panels, then you don't have to exclude improper ancestry.
 
removedWrongAnceInst <- function(plink, inputPrefix, metaDataFile, ancestrySymbol, outputPrefix){

	metaData = read.table(metaDataFile, stringsAsFactors=FALSE, header=TRUE) 
	EAids = metaData[which(metaData[,"ance"]==ancestrySymbol),"IID"]
	fam = read.table(file=paste0(inputPrefix, ".fam"), stringsAsFactors=FALSE) 
	famEA = fam[is.element(fam[,2], EAids), 1:2]
	plinkFormatDat = famEA
	write.table(plinkFormatDat, file=paste0(outputPrefix,".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")  

	system( paste0(plink, " --bfile ", inputPrefix, " --keep ", paste0(outputPrefix,".txt"), " --make-bed --out ", outputPrefix) )
	system( paste0("rm ",  outputPrefix, ".txt") )

}


##########################################   
##########################################
## step 6
## Remove the excluded probes of the chip and check that there are not two probes with the same name afterwards
## > remove AFFX, cnvi, or similar snps
 
removedExclProbe <- function(plink, inputPrefix, excludedProbeIds, outputPrefix){
 
	if (!is.null(excludedProbeIds)) {

		system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", excludedProbeIds, " --make-bed --out ", outputPrefix) )
		## remove .txt
		system( paste0("rm ", excludedProbeIds) )
	} else { 
		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
	    }

}
  
## step7 

 	########################################  outputPrefix >  remove unmapped 2 chipRef
 	## re-generate chip annotation file/prepare chipAnnoFile 
 	## generate chip annotation file in a correct format

removedUnmapProbes <- function(plink, inputPrefix, chipAnnoFile.svn, outputPrefix, outputPrefix.snp){

 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(input=chipAnnoFile.svn, output=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(input=chipAnnoFile.svn, output=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)
	# system( paste0("rm ", annoFile))
	## check the overlapping SNPs
	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
	interSNPs = intersect(bim[,2], chipAnno[,"chipSnpID"])  
	unmapped = setdiff(bim[,2], interSNPs)
	write.table(unmapped, file=paste0(outputPrefix.snp, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix.snp, ".txt"), " --make-bed --out ", outputPrefix) )
 
}



# step8 
## remove double SNPs, dif SNP-A IDs but same rs-names 

 

removedDoubleProbes <- function(plink, inputPrefix, chipAnnoFile.svn, chipType, outputPrefix, outputPrefix.snp){
 	
 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(input=chipAnnoFile.svn, output=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(input=chipAnnoFile.svn, output=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE) 

	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=F) 
	## cbind (combine) bim and chip annotation files  
	chipAnnoV1 = chipAnno[is.element(chipAnno[,"chipSnpID"], bim[,2]),]
	chipAnnoV1sort = chipAnnoV1[match(bim[,2], chipAnnoV1[,"chipSnpID"]),]
	comb = cbind(bim, chipAnnoV1sort)
	############# remove SNPs which have a duplicated rs-name or position (i.e. bp and chr) in this file
	## remove SNPs with duplicated position first
	chrNames = names(table(comb[,1]))
	dupPos = c()
	for (i in chrNames) { 
		print(i)
		subData = comb[which(comb[,1]==i), ]
  		subDup = subData[duplicated(subData[,"pos"]) | duplicated(subData[,"pos"], fromLast=TRUE), ] 
 		dupPos = rbind(dupPos, subDup)
	} 
	# print(dupPos)	
	snpWithdupPos = dupPos[,"chipSnpID"]
	## return all the duplicated rs-names (not only either one)

	if (chipType=="affymetrix") { 
 		snpdup = comb[duplicated(comb[,"rsID"]) | duplicated(comb[,"rsID"], fromLast=TRUE), "chipSnpID"] 
 	} else if (chipType=="illumina"){ 
 		snpdup = comb[duplicated(comb[,"V2"]) | duplicated(comb[,"V2"], fromLast=TRUE), "chipSnpID"]  
 	 }

 	allDupSNPs = c(snpWithdupPos, snpdup) 
 	allDupSNPs = unique(allDupSNPs)
	write.table(allDupSNPs, file=paste0(outputPrefix.snp, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	cmd = paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix.snp, ".txt"), " --make-bed --out ", outputPrefix)  
	system(cmd)  
 
}



## step 9 
###################### updateSNPinfo

updatedSnpInfo <- function(plink, inputPrefix, chipAnnoFile.svn, chipType, inputPrefix.rs, inputPrefix.chr, inputPrefix.pos, inputPrefix.strand, outputPrefix){
 	
 	annoFile = "chipAnnoRefb37.txt"
 	if (chipType=="affymetrix") { 
 		prepareChipAnnoFile4affymetrix(input=chipAnnoFile.svn, output=annoFile)
 	} else if (chipType=="illumina"){ 
 		prepareChipAnnoFile4Illumina(input=chipAnnoFile.svn, output=annoFile)
 	  }
 
	## find the overlapping
	chipAnno = read.table(file=annoFile, header=TRUE, stringsAsFactors=FALSE)
	system( paste0("rm ", annoFile)) ## not used anymore

	# ## cbind (combine) bim and chip annotation files 
	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE) 
	interSNPs = intersect(bim[,2], chipAnno[,"chipSnpID"]) 
	bimV1 = bim[is.element(bim[,2], interSNPs),]
	chipAnnoV1 = chipAnno[is.element(chipAnno[,"chipSnpID"], interSNPs),]
	chipAnnoV1sort = chipAnnoV1[match(bimV1[,2], chipAnnoV1[,"chipSnpID"]),]
	comV2 = cbind(bimV1, chipAnnoV1sort) 
 
	## Update main info  
	if (chipType=="affymetrix") { 
 			updateSNP2rs = subset(comV2, select=c(V2, rsID))
			updateSNPchr = subset(comV2, select=c(rsID, chr))
			updateSNPpos = subset(comV2, select=c(rsID, pos)) 
			updateSNPbackward = comV2[which(comV2[,"strand"]=="-"), "rsID"]  ## note strand 
 	} else if (chipType=="illumina"){ 
	    	updateSNP2rs = subset(comV2, select=c(V2, V2)) ## no need to change but for consistency with Affymtrix
			updateSNPchr = subset(comV2, select=c(V2, chr))
			updateSNPpos = subset(comV2, select=c(V2, pos)) 
			updateSNPbackward = comV2[which(comV2[,"strand"]=="-"), "V2"]  ## note strand 
 	 }

	write.table(updateSNP2rs, file=paste0(inputPrefix.rs, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPchr, file=paste0(inputPrefix.chr, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPpos, file=paste0(inputPrefix.pos, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	write.table(updateSNPbackward, file=paste0(inputPrefix.strand, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	
	## update rs, chr, and pos one by one 	## flip to the forward strand
	system( paste0(plink, " --bfile ", inputPrefix,    " --update-name ", paste0(inputPrefix.rs, ".txt"), " 2 1  --make-bed --out ", inputPrefix.rs) )  
	system( paste0(plink, " --bfile ", inputPrefix.rs,  " --update-chr ", paste0(inputPrefix.chr, ".txt"), " 2 1 --make-bed --out ", inputPrefix.chr) )  
	system( paste0(plink, " --bfile ", inputPrefix.chr, " --update-map ", paste0(inputPrefix.pos, ".txt"), " 2 1 --make-bed --out ", inputPrefix.pos) )   
	system( paste0(plink, " --bfile ", inputPrefix.pos, " --flip ", paste0(inputPrefix.strand, ".txt"), " --make-bed --out ", inputPrefix.strand) )  

	## copy/rename all snp info updated plink files
	system( paste0("cp ", inputPrefix.strand, ".bed ", outputPrefix, ".bed") )
	system( paste0("cp ", inputPrefix.strand, ".bim ", outputPrefix, ".bim") )
	system( paste0("cp ", inputPrefix.strand, ".fam ", outputPrefix, ".fam") )
 
 	## remove all tmp files (rs, chr, pos)
	system( paste0("rm ", inputPrefix.rs, ".*") )  
	system( paste0("rm ", inputPrefix.chr, ".*") )  
	system( paste0("rm ", inputPrefix.pos, ".*") )  
 	system( paste0("rm ", inputPrefix.strand, ".*"))

}


## step 10 

changedXyChr <- function(plink, inputPrefix, outputPrefix){
 	## split X chr into PAR(chr25) and non-PAR (chr23)
 	system( paste0(plink, " --bfile ", inputPrefix, "  --split-x hg19 --make-bed --out ", outputPrefix ))
}

  
## step 11 

removedYMtSnp <- function(plink, inputPrefix, outputPrefix){

 	bim = read.table(paste0(inputPrefix, ".bim"), stringsAsFactors=FALSE)  
 	snpChrY =  bim[which(bim[,1]== 24), 2] ##  
	snpChrMt = bim[which(bim[,1] == 26), 2] 
	snpChrYMt = c(snpChrY, snpChrMt) 
	write.table(snpChrYMt, file=paste0(outputPrefix, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
	system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix) )  
	system( paste0("rm ", paste0(outputPrefix, ".txt")) )
}

  

 
 
#########################################################################################  
######################################################################################### 

prepareChipAnnoFile4affymetrix <- function(input, output){

	inputNew = paste0(input, "New")
	system( paste0("sed 1d ", input, " > ", inputNew) )
	chipAnnoRefraw = read.table(file=inputNew, stringsAsFactors=FALSE) 
	# colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "strand") ## for Illumina
	colnames(chipAnnoRefraw) = c("chipSnpID", "rsID", "chr", "pos", "strand")

	## remove SNPs with strange strand 
	whUnknown = which(chipAnnoRefraw[,"strand"]=="---") 
	chipAnnoRefraw2 = chipAnnoRefraw[-whUnknown,]

	## only see 3 different cases (if 25--> XY)
	whX = which(chipAnnoRefraw2[,"chr"]=="X")
	whY = which(chipAnnoRefraw2[,"chr"]=="Y")
	whMT = which(chipAnnoRefraw2[,"chr"]=="MT")

	chipAnnoRefraw2[whX,"chr"] = 23 
	chipAnnoRefraw2[whY,"chr"] = 24
	chipAnnoRefraw2[whMT,"chr"] = 26
  
  	whAFF = grep("AFFX-SNP", chipAnnoRefraw2[,"chipSnpID"])
	# str(whAFF)
	chipAnnoRefraw3 = chipAnnoRefraw2[-whAFF,]
	write.table(chipAnnoRefraw3, file=output, quote=F, row.names=F, col.names=TRUE, eol="\r\n", sep="\t")

}


#########################################################################################  
######################################################################################### 
## prepareChipAnnoFile4Illumina

# colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "match2genom", "strand", "allele")
## the column name must be fixed
## chipSnpID is set in the way so that it's comparable to Affymetrix chip

## 1. don't have to remove "cnvi", since we will remove this type of SNPs in .raw plink files
## 2. don't have to remove 'unknown' SNPs, since they won't map to our raw plink files.
 
prepareChipAnnoFile4Illumina <- function(input, output){

	chipAnnoRefraw = read.table(file=input, stringsAsFactors=FALSE)
	# print(colnames(chipAnnoRefraw))
	## c("chipSnpID", "chr", "pos", "match2genom", "strand", "allele")
	chipAnnoRefraw <- chipAnnoRefraw[,-c(4, 6)]
	colnames(chipAnnoRefraw) = c("chipSnpID", "chr", "pos", "strand")

	## only see 3 different cases (if 25--> XY)
	whX = which(chipAnnoRefraw[,"chr"]=="X")
	whY = which(chipAnnoRefraw[,"chr"]=="Y")
	whMT = which(chipAnnoRefraw[,"chr"]=="MT")

	chipAnnoRefraw[whX,"chr"] = 23 
	chipAnnoRefraw[whY,"chr"] = 24
	chipAnnoRefraw[whMT,"chr"] = 26 

	chipAnnoNew = chipAnnoRefraw
	write.table(chipAnnoNew, file=output, quote=F, row.names=F, col.names=TRUE, eol="\r\n", sep="\t")

}

  