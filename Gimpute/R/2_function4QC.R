 
  
# File   : 2_function4QC.R
# Author : Junfang Chen
# Version0: 22 Jul 2016
# VersionX: 08 Jan 2018
  
## This script is used to define a specific QC pipeline, which is adjusted from standard QC cited from 108 loci nature paper.
## Population strucutre is detected via PCA. The cutoff for outlier detection is observed by looking at the generated PCA plot. 
  

 
##########################################   
##########################################

# Description:
# Determine for each SNP of the chromosome 23 from the genotype data
# the number of male instances which have the value one as the minor
# allele count for that SNP and remove all SNPs which number is higher
# than 0.5 % of the number of male instances.


# Arguments:
# outputPrefix: Primary output files:  
# Output file with the number of instances with heterozygous alleles for
# each SNP of the chromosome 23 before SNP removal (each line contains
# a SNP name and the respective number, lines are sorted descending by
# number):
# outputPrefix.snp: 
# Output file with the number of instances with heterozygous alleles for
# each SNP of the chromosome 23 after SNP removal:


removedSnpHetX <- function(plink, inputPrefix, cutoff, outputPrefix, outputPrefix.snp){

		## in order to remove heterozygous SNPs for male in chrX, two extra steps are required.
		## 2_01_removedHetXSnp.%

		## just to get .hh file and .fam file 
		system( paste0(plink, " --bfile ", inputPrefix, " --chr 23 --filter-males --make-bed --out male23nonPAR" ))

		hh = read.table("male23nonPAR.hh", stringsAsFactors=F)
		fam = read.table("male23nonPAR.fam", stringsAsFactors=F)

		hetSNPsFreq = table(hh[,3])
		# hetSNPFreqFreq = table(hetSNPs)

		cutoff4removeHetSNP = nrow(fam)*cutoff
		mostFakeSNPs = hetSNPsFreq[which(hetSNPsFreq>= cutoff4removeHetSNP)] 
		mostFakeSNPs = names(mostFakeSNPs)
		write.table(mostFakeSNPs, file="mostFakeSNPs.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
		## remove these fake SNPs
		system( paste0(plink, " --bfile ", inputPrefix, " --exclude mostFakeSNPs.txt --make-bed --out ", outputPrefix) )
		system( paste0("rm mostFakeSNPs.txt") )
		system( paste0("rm male23nonPAR.*") )
		system( paste0("rm ", inputPrefix,".*") )

		## generate hetSNPsFreq in .txt file 
		hetSNPsWithInstNumber = data.frame(hetSNPsFreq, stringsAsFactors=F)
		hetSNPsWithInstNumber = hetSNPsWithInstNumber[order(hetSNPsWithInstNumber[,2], decreasing=T),] 
		write.table(hetSNPsWithInstNumber, file=paste0(outputPrefix.snp, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
}

 
 
##########################################   
##########################################



removedInstHetX <- function(plink, inputPrefix, cutoff, outputPrefix, outputPrefix.snp){

	if ( is.na(file.size(paste0(inputPrefix, ".hh")))==TRUE ){ 

		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
		system( paste0("touch ", outputPrefix.snp, ".txt") ) 
  
	} else { 	

		hh = read.table(paste0(inputPrefix, ".hh"), stringsAsFactors=F)
		fam = read.table(paste0(inputPrefix, ".fam"), stringsAsFactors=F)

		hetInstFreq = table(hh[,2])
		# table(hetInstFreq) 
		cutoff4removeHetInst = cutoff
		mostFakeInst = hetInstFreq[which(hetInstFreq>= cutoff4removeHetInst)] 
		mostFakeInst4plink = fam[is.element(fam[,2], names(mostFakeInst)), 1:2]

		write.table(mostFakeInst4plink, " --bfile ", file="mostFakeInst4plink.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
		## remove these fake SNPs
		system( paste0(plink, " --bfile ", inputPrefix, " --remove mostFakeInst4plink.txt --make-bed --out ", outputPrefix) )
		system( paste0("rm mostFakeInst4plink.txt") )

		## generate hetSNPsFreq in .txt file 
		InstWithHetSNPs = data.frame(hetInstFreq, stringsAsFactors=F)
		## sort -nr -k 2 
		InstWithHetSNPs = InstWithHetSNPs[order(InstWithHetSNPs[,2], decreasing=T),] 
		write.table(InstWithHetSNPs, file=paste0(outputPrefix.snp, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 

	}	
}




##########################################   
##########################################
## Set all heterozygous alleles of SNPs of the chromosome 23 for males
# (i.e. when the SNP has the value one as the minor allele count) as
# missing.



setHeteroHaploMissing <- function(plink, inputPrefix, outputPrefix){

	# system( paste0(plink, " --bfile ", inputPrefix, " --set-hh-missing --make-bed --out 3tmp") )  
	system( paste0(plink, " --bfile ", inputPrefix, " --set-hh-missing --make-bed --out ", outputPrefix) )  

	# ## copy/rename all snp info updated plink files
	# system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
	# system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
	# system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
}  




 


##########################################   
##########################################

removedSnpMissPre <- function(plink, snpMissing, inputPrefix, outputPrefix){
	if (!is.null(snpMissing)) { 
		system( paste0(plink, " --bfile ", inputPrefix, " --geno ", snpMissing, " --make-bed --out ", outputPrefix) )  
	}	
}  
  

##########################################   
########################################## 
# subject missingness < 0.02; 
removedInstMiss <- function(plink, sampleMissing, inputPrefix, outputPrefix){
	if (!is.null(sampleMissing)) {
		system( paste0(plink, " --bfile ", inputPrefix, " --mind ", sampleMissing, " --make-bed --out ", outputPrefix) ) 
		system( paste0("rm ", outputPrefix, ".irem") )
 	}
}   


##########################################   
##########################################	

# 6. autosomal heterozygosity deviation (| Fhet | < 0.2); 
##  This analysis will automatically skip haploid markers (male X and Y chromosome markers).

removedInstFhet <- function(plink, Fhet, inputPrefix, outputPrefix){ 

	if (!is.null(Fhet)) {
		system( paste0(plink, " --bfile ", inputPrefix, " --het --out ", outputPrefix) )
		#  F inbreeding coefficient estimate
		autoHet = read.table(file=paste0(outputPrefix, ".het"), header=T)  
		fhet = autoHet[, "F"]
		qc_data_fhet = autoHet[which(abs(fhet) >= Fhet), c(1, 2)]  
		## the individual IDs to be removed  
		write.table(qc_data_fhet, file=paste0(outputPrefix, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		##############
		# To remove certain individuals 
		system( paste0(plink, " --bfile ", inputPrefix, " --remove ", paste0(outputPrefix, ".txt"), " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", outputPrefix, ".het") )
		system( paste0("rm ", outputPrefix, ".txt") )
 	} 
}   

##########################################   
##########################################

## plink2 --bfile unpruned_data --make-king-table --king-cutoff 0.177 --king-table-filter 0.177 --make-bed --out pruned_data
## https://groups.google.com/forum/#!topic/plink2-users/F-b4XRF8CSc
  
removedInstRelated <- function(plink, kinshipValue, inputPrefix, outputPrefix, outputPrefix.ID){ 

	# new step6. 6. Replace the paternal ID and maternal ID of instances  by the value zero  
	### special case for plink 
	# plink2 = '/data/noether/tools/plink/plink2'
	# system( paste0(plink2, " --bfile ", inputPrefix, " --make-king-table --king-cutoff 0.177  --king-table-filter 0.11 --make-bed --out ", outputPrefix) ) 
	## >> this doesn't work; 

	## copy/rename all snp info updated plink files
	system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
	system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
	system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
	system( paste0("touch ", outputPrefix.ID, ".txt") )  ## empty

}   
	 

 

##########################################   
##########################################
 

# Replace the paternal ID and maternal ID of instances (childs) by the
# value zero if the paternal ID and the maternal ID do not belong to any
# instance (parent) with the same family ID as the child. 


removedParentIdsMiss <- function(inputPrefix, outputPrefix){ 

	# Remove the parent IDs which do not belong to instances.
	system( paste0(plink, " --bfile ", inputPrefix, " --make-founders require-2-missing --make-bed --out ", outputPrefix) ) 

	# 	## copy/rename all snp info updated plink files
	# system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
	# system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
	# system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )
}   



##########################################   
##########################################  should also include case/control

 # SNP missingness < 0.02 (after sample removal); >> only for case and control data >> so keep as step 9
 
removedSnpMissPost <- function(plink, snpMissing2, inputPrefix, outputPrefix){ 

	if (!is.null(snpMissing2)) {
		system( paste0(plink, " --bfile ", inputPrefix, " --geno ", snpMissing2, " --make-bed --out ", outputPrefix) )  
	} 

}   

  # testMissing = NULL # 7  !!! ! ## different from case control data  >>>>> tobeImproved


##########################################   
##########################################  

# SNP missingness < 0.02 (after sample removal);
removedSnpMissDiff <- function(plink, snpMissing2, inputPrefix, outputPrefix){ 

	system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
	system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
	system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

}   



##########################################   
##########################################  
#.Remove chrX SNPs with (missingness >= 0.05) in females
 
removedSnpFemaleChrXmiss <- function(plink, femaleChrXmiss, inputPrefix, outputPrefix){ 
 
 	if ( femaleChrXmiss==TRUE ) {
 	### additional QC (female-chrX SNPs, missingness ok?)  
 	 	outputPrefix.tmp1 = paste0(outputPrefix, "tmp1")
		outputPrefix.tmp2 = paste0(outputPrefix, "tmp2")
		system( paste0(plink, " --bfile ", inputPrefix, " --filter-females --chr 23 --make-bed --out ", outputPrefix.tmp1) )
		system( paste0(plink, " --bfile ", inputPrefix, " --filter-females --chr 23 --geno 0.05 --make-bed --out ", outputPrefix.tmp2) )
 	 
		 ## check if equal  
		femaleChrXorig = read.table(paste0(outputPrefix.tmp1, ".bim"), stringsAsFactors=FALSE) 
		femaleChrXMiss = read.table(paste0(outputPrefix.tmp2, ".bim"), stringsAsFactors=FALSE)  
		snps2removed = setdiff(femaleChrXorig[,2], femaleChrXMiss[,2]) 
		snps2removedfile = paste0(outputPrefix, ".txt")
		write.table(snps2removed, file=snps2removedfile, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 

		system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", snps2removedfile, " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", outputPrefix.tmp1, ".*") )
		system( paste0("rm ", outputPrefix.tmp2, ".*") )
		system( paste0("rm ", snps2removedfile) ) 
    }
}


##########################################   
##########################################

## Remove autosomal SNPs with HWE (p < 10?6) in controls.
removedSnpHWEautoControl <- function(plink, hweAUTOcontrol, pval, inputPrefix, outputPrefix){ 
	   
	#  . Remove autosomal SNPs with HWE (p < 10−6) in controls.
	if ( hweAUTOcontrol==TRUE ) {

		outputPrefix.tmp = paste0(outputPrefix, "tmp")
		system(  paste0(plink, " --bfile ", inputPrefix, " --filter-controls --hardy --autosome --make-bed --out ", outputPrefix.tmp) )  
		###  read HWE p values 
		hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), header=T, stringsAsFactors=F) 
		hweControl = hweCheck[which(hweCheck$TEST == "UNAFF"), ] # ## for controls 
		snpHweValuesAutoCt = subset(hweControl, select=c(SNP, P))
		write.table(snpHweValuesAutoCt, file=paste0(outputPrefix.pval, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		removedSNPs_hweControl = hweControl[which(hweControl$P <= pval), "SNP"]
		write.table(removedSNPs_hweControl, file=paste0(outputPrefix.snp, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

		########## exclude SNPs 
		system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix.snp, ".txt"), " --make-bed --out ", outputPrefix) ) 
		system( paste0("rm ", outputPrefix.tmp, ".*") )
 	}

}


##########################################   
##########################################  
## Remove chrX SNPs with HWE (p < 10.6) in female controls. 
removedSnpFemaleChrXhweControl <- function(plink, femaleChrXhweControl, pval, inputPrefix, outputPrefix){ 
	    
 	#  . Remove chrX SNPs with HWE (p < 10−6) in female controls.
 	if ( femaleChrXhweControl==TRUE ) {
  
		outputPrefix.tmp = paste0(outputPrefix, "tmp") 
		system( paste0(plink, " --bfile ", inputPrefix, " --filter-females --filter-controls --chr 23 --hardy ", " --make-bed --out ", outputPrefix.tmp) )
		## read p values
		hweCheck <- read.table(file=paste0(outputPrefix.tmp, ".hwe"), header=T, stringsAsFactors=F) 
		hweControl = hweCheck[which(hweCheck$TEST == "UNAFF"), ] # ## for controls 
		snpHweValuesChrXCt = subset(hweControl, select=c(SNP, P))
		write.table(snpHweValuesChrXCt, file=paste0(outputPrefix.pval, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		## remove bad SNPs
		removedSNPs_hweControl = hweControl[which(hweControl$P <= 0.000001), "SNP"]
		write.table(removedSNPs_hweControl, file=paste0(outputPrefix.snp, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")

		system( paste0(plink, " --bfile ", inputPrefix, " --exclude ", paste0(outputPrefix.snp, ".txt"), " --make-bed --out ", outputPrefix) )
		system( paste0("rm ", outputPrefix.tmp, ".*") ) 
	}
 
}
 

## step 14
## for genotyping data with only control IDs, step 14 is ingored. >> or it is the same as step 13 
## all female instances == all female controls

## step 15
## for genotyping data with only control IDs, step 14 is ingored. >> or it is the same as step 13 
## no case IDs
	  

## http://cnsgenomics.com/software/gcta/#Download
##########################################   
##########################################
## required libraries  
library(lattice) ## PCA plot

plotPCA4plink <- function(gcta, inputPrefix, outputPCs){ 

	autosomefn = paste0(inputPrefix, "Autosome")
	system( paste0(gcta, " --bfile ", inputPrefix, " --make-grm --autosome --out ", autosomefn, " --thread-num 30") )
	system( paste0(gcta, " --grm ", autosomefn, " --pca 20 --out ", autosomefn, " --thread-num 30") )

	eigen = read.table(file=paste0(autosomefn,".eigenvec"), stringsAsFactors=F)
	pcs = eigen[,1:4] ## first two PCs
	write.table(pcs, paste0(outputPCs, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")
	pcWithGroup = cbind(pcs, stringsAsFactors=F)

	png(paste0(outputPCs, ".png"), width=8, height=6, units="in", res=800)
	print( xyplot(pcWithGroup[,4] ~ pcWithGroup[,3], data=pcWithGroup, 
	       auto.key=list(space="right"),  
	       jitter.x=TRUE, jitter.y=TRUE, xlab="PC1", ylab="PC2") )
	dev.off()
	## remove unwanted files
	system( paste0("rm ", autosomefn, ".*") )
}


 
##
######################################################
######################################################
 
## remove selfAA  >> is removed in the first section 1-conversion
 
removeOutlierByPCs <- function(inputPCs, cutoff, cutoffSign, outputIDs, plink, gcta, inputPrefix, outputPCplot, outputPrefix) {

	 
	## if no outliers or no need to remove PC outliers. 
	if ( is.null(cutoff) ==TRUE ){ 

		## copy/rename all snp info updated plink files
		system( paste0("cp ", inputPrefix, ".bed ", outputPrefix, ".bed") )
		system( paste0("cp ", inputPrefix, ".bim ", outputPrefix, ".bim") )
		system( paste0("cp ", inputPrefix, ".fam ", outputPrefix, ".fam") )

	} 
			
	else {  

		PCs =read.table(file=paste0(inputPCs, ".txt"), stringsAsFactors=FALSE)   
		if (length(cutoff) > 1) { ## if the outliers should be removed on both side of the cluster
				outliersPC1v1 = PCs[which(PCs[,3] <= cutoff[1]), 1:2] ## detected by PC1
				outliersPC1v2 = PCs[which(PCs[,3] >= cutoff[2]), 1:2] ## detected by PC1
				outliersPC1 = rbind(outliersPC1v1, outliersPC1v2)

			} else {
			  	  if (cutoffSign == "smaller"){ 
			  	  	outliersPC1 = PCs[which(PCs[,3] <= cutoff), 1:2] ## detected by PC1
			  	  	} else if (cutoffSign == "greater"){
			  	  		outliersPC1 = PCs[which(PCs[,3] >= cutoff), 1:2] ## detected by PC1
			  	  	 }
				}	 
	 
		write.table(outliersPC1, paste0(outputIDs, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")

		system( paste0(plink, " --bfile ", inputPrefix, " --remove ", outputIDs, ".txt --make-bed --out ", outputPrefix) )

		## Plot first two PCs again
		plotPCA4plink(gcta, outputPrefix, outputPCplot)

		outliersPC1sampleID = outliersPC1[,2] ## re-write the output --> so only only sample ID generated
		write.table(outliersPC1sampleID, paste0(outputIDs, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")

	}	
}
 

## tobeImproved for the last step in the above
# Output file with the IDs of the removed instances and the values of
# their eigenvectors, sorted ascending by the PC values:







