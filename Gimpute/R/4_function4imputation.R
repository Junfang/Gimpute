
# File   : 4_function4imptation.R
# Author : Junfang Chen
# Version0: 06 Jun 2016
# VersionX: 15 Jan 2018
 

##########################################################################
########################################################################## removedMonoSnp.R
# Description:
# Remove monomorphic SNPs(double check if mismatched alleles) from the Lifted and QC instance data.
# Arguments:
# inputPrefix: Prefix of the input plink files (after QC and lifting)
# inputPrefix: Prefix of the output plink files (after removing monomorphic SNPs)


removedMonoSnp <- function(plink, inputPrefix, outputPrefix, outputPrefix.snp, outputPrefixTmp){  

		## copy plink files
		tmpFile = "gwas_dataQCtmp"
		system( paste0("cp ./3-lifting/", inputPrefix, ".bed  ./4-imputation/tmp4impute/1-checkAlign/", tmpFile, ".bed") )
		system( paste0("cp ./3-lifting/", inputPrefix, ".fam  ./4-imputation/tmp4impute/1-checkAlign/", tmpFile, ".fam") )
		system( paste0("cp ./3-lifting/", inputPrefix, ".bim  ./4-imputation/tmp4impute/1-checkAlign/", tmpFile, ".bim") )
		  
		## go to   
		setwd("./4-imputation")
		setwd("./tmp4impute" )  
		setwd("./1-checkAlign")
 
		## input  
		bim = read.table(paste0(tmpFile, ".bim"), stringsAsFactors=F)
		monoSNPs = bim[which(bim[,5]==0),2] 
		outputPrefix.snpTmp = "monoSNPs.txt"
		write.table(monoSNPs, file=outputPrefix.snpTmp, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ") 
		setwd("..")
		 
		system(paste0("mv ./1-checkAlign/", tmpFile, ".fam ./2-dataFiles/"))
		system(paste0("mv ./1-checkAlign/", tmpFile, ".bed ./2-dataFiles/"))
		system(paste0("mv ./1-checkAlign/", tmpFile, ".bim ./2-dataFiles/"))
		system(paste0("mv ./1-checkAlign/", outputPrefix.snpTmp, " ./2-dataFiles/"))
		 
		setwd("./2-dataFiles/") 
		## exclude those not well aligned snps  ## intermediate output 
		system(paste0(plink, " --bfile ", tmpFile, " --exclude ", outputPrefix.snpTmp, " --make-bed --out ", outputPrefixTmp)) 
		# system( paste0("rm *.log *.hh")) ## *.hh may not exist any more

		setwd("..") ##   /tmp4impute
		setwd("..") ##   /4-imputation
		
		# copy final removed plink files and the file with removed SNPs
		system( paste0("scp ./tmp4impute/2-dataFiles/", outputPrefixTmp, ".bed ", outputPrefix, ".bed"))
		system( paste0("scp ./tmp4impute/2-dataFiles/", outputPrefixTmp, ".fam ", outputPrefix, ".fam"))
		system( paste0("scp ./tmp4impute/2-dataFiles/", outputPrefixTmp, ".bim ", outputPrefix, ".bim"))
		system( paste0("scp ./tmp4impute/2-dataFiles/", outputPrefix.snpTmp, " ", outputPrefix.snp, ".txt"))

}





##########################################################################
########################################################################## Main imputation pipeline
 
# Use the imputation reference files to generate pre-phase haplotypes by
# SHAPEIT and then do the imputation by using IMPUTE2.

## sub-steps
# step 2.1 chrWiseSplit.R
# step 2.2 chunk4eachChr.R
# step 2.3 prePhasingByShapeit.R
# step 2.4 imputedByImpute2.R 
# step 2.5 formatConvertGtool.R 
# step 2.6 mergeImputeData.R 
# step 2.7 filterImputeData.R 


# >>>

##########################################################################
########################################################################## chrWiseSplit.R
# Description:
# split the whole genome genotyping data chromosome-wise; allow parallel computating for all chromosomes.
# if chromosome 25 is also available, further split chr25 (PAR or Chr_XY) into PAR1 and PAR2 according to the genomic coordination GRCh37 

# from https://en.wikipedia.org/wiki/Pseudoautosomal_region.
# The locations of the PARs within GRCh37 are:  
# PAR1	X	60001	2699520 
# PAR2	X	154931044	155260560 

# Arguments:
# inputPrefix: Prefix of the input plink files  (before splitting).
# 			   Also it's the Prefix of the output plink files (chr-wise genotyping plink files).
# X_PAR1sufix:  if chr 25 is available and with PAR1, then generate sufix with X_PAR1 for chrX_PAR1 
# X_PAR2sufix:  if chr 25 is available and with PAR2, then generate sufix with X_PAR2 for chrX_PAR2 


library(doParallel) 
chrWiseSplit <- function(plink, inputPrefix, X_PAR1sufix, X_PAR2sufix){ 

	setwd("./tmp4impute") ## 
	setwd("./2-dataFiles/")
   
    ## check which chromosomes are available to be splitted from the .bim file
	bim = read.table(paste0(inputPrefix, '.bim'), stringsAsFactors=F)
	chrs = as.integer(names(table(bim[,1])))
 
	chrslist = as.list(chrs)
	mclapply(chrslist, function(i){
		cmd = paste0(plink, " --bfile ", inputPrefix, " --chr ", i, " --make-bed --out ", inputPrefix, i)  
		system(cmd)
	}, mc.cores=length(chrs))

	## if chromosome 25 is also available then re-arrange it
 	if (is.element(25, chrs)){  

 		print('PAR is available in chrX!') 
		bim25 = read.table(paste0(inputPrefix, "25.bim"), stringsAsFactors=F) 
		pos4PAR1= c(60001, 2699520) 
		## first check for PAR1 and afterwards for PAR2
		if ( length(which(bim25[,4]<=pos4PAR1[2]))!=0 ){ 
	   
	   		print('PAR1 is available in chrX!')
	   		bimPos4par1 = which(bim25[,4]<=pos4PAR1[2])
			rs4PAR1 = bim25[bimPos4par1,2]
			write.table(rs4PAR1, file="rs4PAR1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
			system( paste0(plink, " --bfile ", inputPrefix, "25 --extract rs4PAR1.txt --make-bed --out ", inputPrefix, X_PAR1sufix) )
			par1 = TRUE  
			## check for PAR2, if any SNPs out of PAR1, then PAR2 also available 
			if (length(bimPos4par1)<nrow(bim25)){
				print('PAR2 is available in chrX!') 
				## just to exclude rs4PAR1.txt
				system( paste0(plink, " --bfile ", inputPrefix, "25 --exclude rs4PAR1.txt --make-bed --out ", inputPrefix, X_PAR2sufix) )
				par2 = TRUE 
			} else {print('PAR2 is NOT available in chrX!') }

		} else { 
			print('PAR2 is available in chrX! But NOT PAR1, all chr25 on PAR2')
			par1 = FALSE
			par2 = TRUE
			system( paste0("cp ", inputPrefix, "25.bed ", inputPrefix, X_PAR2sufix, ".bed") )
			system( paste0("cp ", inputPrefix, "25.bim ", inputPrefix, X_PAR2sufix, ".bim") )
			system( paste0("cp ", inputPrefix, "25.fam ", inputPrefix, X_PAR2sufix, ".fam") )
		} 
	} else {  
		print('PAR is NOT available in chrX!') 
		par1 = FALSE
		par2 = FALSE
	} 

	return(par=list(par1, par2))
	## go to the main directory > ./tmp4impute
	setwd("..")   
}
 


##########################################################################
########################################################################## chunk4eachChr.R
# Description:
# chunking each chromosome into multiple segments by a predefined window size

# Arguments:
# inputPrefix: Prefix of the input plink files for each chromosome.
# outputPrefix: Prefix of the output plink files for each chunk 
# chrs: specifiy the chromosomes for chunking
# windowSize: the window size for each chunk   
  

chunk4eachChr <- function(inputPrefix, outputPrefix, chrs, windowSize){  

		for (i in chrs){ 

			bimfilename = paste0("./2-dataFiles/", inputPrefix, i, ".bim")
			bimdata = read.table(file=bimfilename, sep="\t", stringsAsFactors=F)
			position = bimdata[,4]
			posStart = head(position,1)
			posEnd = tail(position,1)
			chunkStart = seq(posStart, posEnd, windowSize)
			chunkEnd = chunkStart + windowSize -1
			chunkMatrix = cbind(chunkStart, chunkEnd)

			## positions are only within a chunk
			if (nrow(chunkMatrix) == 1){
				chunks = chunkMatrix
			} else {  
				##  it may happen that only a few SNPs from the last chunk; but if the last chunk is large, then specify -allow_large_regions 
				chunks = head(chunkMatrix, -1) ## merge last-second to last  
				chunks[nrow(chunks), 2] <- posEnd
			}
			

			## check if any chunk with NO snps within it   
			SNPcountsPerChunk <- c() 
			for (j in 1:nrow(chunks)){
				chunkbottom = chunks[j,1]
				chunkup = chunks[j,2]
				tmp = which(position >= chunkbottom & position <= chunkup)  ## which fall within chunk
				tmp = length(tmp)
				SNPcountsPerChunk <- c(SNPcountsPerChunk, tmp) 
			} 


		 	wh0 = which(SNPcountsPerChunk==0)
		 	print(paste0("chr",i))
		 	print(wh0)
		 	chunkLength = nrow(chunks) - length(wh0)
		 	## remove such chunks if wh0  

		 	if (length(wh0)!=0){ chunks = chunks[-wh0,]  }
		 	print(nrow(chunks) == chunkLength)

			chunkfilename = paste0("./4-chunkFile/", outputPrefix, i, ".txt")
			write.table(chunks, file=chunkfilename, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		}
}




	

##########################################################################
########################################################################## prePhasingByShapeit.R
# Description:
# chunking each chromosome into multiple segments by a predefined window size

# Arguments:
# chrs: specifiy the chromosomes for phasing.
# nThreat: the number of threads
# effectiveSize: population effective size 
# chrXflag: --chrX flag, specifically for chrX imputation

library(doParallel)
prePhasingByShapeit <- function(prePhasingScript, chrs, nThreat, effectiveSize, chrXflag){

	## prePhasing for the autosome (include PAR)
	## Exclude chr23 at first	
	autoChrs = setdiff(chrs, '23')
	chrslist = as.list(autoChrs)
	mclapply(chrslist, function(i){
		cmd = paste0(prePhasingScript, " ", i, " ", nThreat, " ", effectiveSize, " ", chrXflag=NULL) 
		system(cmd)
	}, mc.cores=1)	

	if (is.element('23', chrs)){ 
		## input the script: prePhasing4chrX  
		system( paste0(prePhasingScript, " 23 ",   nThreat, " ", effectiveSize, " ", chrXflag) ) 
	}  
} 

 
  

##########################################################################
########################################################################## imputedByImpute2.R
# Description:
 
# chrX 
# http://www.shapeit.fr/pages/m03_phasing/imputation.html
# Two comments:
# 1. You must use the -chrX flag for IMPUTE2 to proceed with X chromosome imputation
# 2. You must give the SAMPLE file generated by SHAPEIT to IMPUTE2. 
# This SAMPLE has a sex column that gives the gender of the GWAS individuals.
 


# Arguments:
# chrs: specifiy the chromosomes for phasing.
# nThreat: the number of threads
# effectiveSize: population effective size 
# XPAR: --chrX flag, specifically for chrX imputation


 
# ## for

library(doParallel)
imputedByImpute2 <- function(impute2script, chrs, effectiveSize, nCore, XPAR){

	## prePhasing for the autosome (include PAR)
	## Exclude chr23 at first	
	pureAutoChrs = setdiff(chrs, c("23", "X_PAR1", "X_PAR2")) 
	for (i in pureAutoChrs){ 	
		chunkfn = paste0("./4-chunkFile/chunks_chr", i, ".txt")
		chunks = read.table(chunkfn, sep=" ")
		dim(chunks) 
		 
		chunklist = as.list(1:nrow(chunks))
		mclapply(chunklist, function(j){
			## imputation scripts
		 	systemCall <- paste0(impute2script, " ", i, " ", chunks[j,1], " ", chunks[j,2], " ", effectiveSize, " ", XPAR=NULL)
			system(systemCall)
		}, mc.cores=nCore)  ## for 64G: it's recommended to use 15 or less cores.  
	}	


	if (is.element(par1, chrs)){ 
		chunkfn = paste0("./4-chunkFile/chunks_chrX_PAR1.txt")
		chunks = read.table(chunkfn, sep=" ")
		dim(chunks)  
		chunklist = as.list(1:nrow(chunks))
		mclapply(chunklist, function(j){
			system( paste0(impute2script, " X_PAR1 ", chunks[j,1], " ", chunks[j,2], " ", effectiveSize, " ", XPAR) )  
		}, mc.cores=nCore) 

	}  

	if (is.element(par2, chrs)){ 
		chunkfn = paste0("./4-chunkFile/chunks_chrX_PAR2.txt")
		chunks = read.table(chunkfn, sep=" ")
		dim(chunks)  
		chunklist = as.list(1:nrow(chunks))
		mclapply(chunklist, function(j){
			system( paste0(impute2script, " X_PAR2 ", chunks[j,1], " ", chunks[j,2], " ", effectiveSize, " ", XPAR) ) 
		}, mc.cores=nCore) 	
	}  


	if (is.element('23', chrs)){  
		chunkfn = paste0("./4-chunkFile/chunks_chr23.txt")
		chunks = read.table(chunkfn, sep=" ")
		dim(chunks)  
		chunklist = as.list(1:nrow(chunks))
		mclapply(chunklist, function(j){
			system( paste0(impute2script4chrX, " 23 ", chunks[j,1], " ", chunks[j,2], " ", effectiveSize) ) 
		}, mc.cores=nCore) 
	}  
} 

 


 

##########################################################################
########################################################################## formatConvertGtool.R
# Description:

# First only extract 'rs' variants after imputation; 
# Convert all chunks of IMPUTE2 format into GTOOL format
 
# Arguments:
# variantPrefix: specify which variants should be retained; in this case the variants start with 'rs'. Namely no INDELs.
# gtoolScript: use Gtool to convert IMPUTE2 format files into PLINK format
# chrs: specifiy the chromosomes for conversion.  
# nCore: the number of cores used for computation.

  
library(doParallel)
formatConvertGtool <- function(variantPrefix, chrs, gtoolScript, nCore){
 		 
		## extract only SNPs starting with 'rs' 
		setwd("./5-imputeResults")  
		ls = system("ls gwas*.impute2", intern=T)
		library(doParallel)  
		biglists = as.list(ls)
		mclapply(biglists, function(i){
			arg1 = paste0(i, "noINDEL.impute2")
			# arg2 = paste0("grep 'rs' ", i, " | awk '{if(length($4)==1 && length($5)==1) print}' | awk '{ if (a[$2]++ == 0) print $0; }' > ", arg1)
			arg2 = paste0("grep '", variantPrefix, "' ", i, " | awk '{if(length($4)==1 && length($5)==1) print}' > ", arg1)
			print(arg2)
			system(arg2)
		}, mc.cores= 38) 
		  

		setwd("..")  
		for (i in chrs){ 
			chunkfn = paste0("./4-chunkFile/chunks_chr", i, ".txt")
			chunks = read.table(chunkfn, sep=" ")
			dim(chunks) 
			chunklist = as.list(1:nrow(chunks))
			mclapply(chunklist, function(j){
			 	cmd <- paste0(gtoolScript, " ", i, " ", chunks[j,1], " ", chunks[j,2])
				system(cmd)
			}, mc.cores=nCore)
		} 
} 






##########################################################################
########################################################################## mergeImputeData.R
# Description:

# merge many chunk.ped chunk.map files together into chr-wise ped and map files!
# create a file containing a list chunk ped and map file names
# then combine all chrs (combine the first 23 chrs; then Xpar)

# Arguments:
 

library(doParallel) 
mergeImputeData <- function(plink, chrs){ 


	setwd('./6-postImpute/')
	## firstly, only consider chromosomes from 1:23; as Xpar chrs are slightly different for processing.
	pureAutoChrs = setdiff(chrs, c("X_PAR1", "X_PAR2")) 
	chrslist = as.list(pureAutoChrs)   
	mclapply(chrslist, function(i){
		pedFile_chr = system(paste0("ls gwas_data_chr", i, ".*.ped"), intern=TRUE)
		mapFile_chr = system(paste0("ls gwas_data_chr", i, ".*.map"), intern=TRUE)	
		pedmap_chr = paste0(pedFile_chr, " ", mapFile_chr)
		fA = gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged = pedmap_chr[-1]
		filesetname = paste0("fileset_chr", i, ".txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		fnew = paste0("gwasImputed_chr", i)

	    arg <- paste0(plink, " --file ", fA, " --merge-list ", filesetname, " --make-bed --out ", fnew)
	    system(arg) 
	}, mc.cores=25)


	############################################################################## combine chrX_PAR and convert into chr25
	if (is.element(c("X_PAR1", "X_PAR2"), chrs)){  
		pedFile_chr = system(paste0("ls gwas_data_chrX_PAR*.ped"), intern=TRUE)
		mapFile_chr = system(paste0("ls gwas_data_chrX_PAR*.map"), intern=TRUE)	
		pedmap_chr = paste0(pedFile_chr, " ", mapFile_chr)
		fA = gsub(".ped", "", pedFile_chr[1])
		pedmap_tobeMerged = pedmap_chr[-1]
		filesetname = paste0("fileset_chr25.txt")
		write.table(pedmap_tobeMerged, file=filesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
		 
		arg <- paste0(plink, " --file ", fA, " --merge-list ", filesetname, " --allow-extra-chr --make-bed --out gwasImputed_oldchr25")
		system(arg) 
		## update chr code for XPAR --> 25
		bim = read.table("gwasImputed_oldchr25.bim", stringsAsFactors=F)
		updateSNPchr = cbind(bim[,2], rep(25, length=nrow(bim))) 
		write.table(updateSNPchr, file="gwasImputed_newchr25.txt", quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ") 
		system( paste0(plink, " --bfile gwasImputed_oldchr25 --allow-extra-chr --update-chr gwasImputed_newchr25.txt 2 1 --make-bed --out gwasImputed_chr25") )  

	}	 
	

	############################################################################## combine all bed files
	bedFile_chr = system(paste0("ls gwasImputed_chr*.bed"), intern=TRUE)
	bimFile_chr = system(paste0("ls gwasImputed_chr*.bim"), intern=TRUE)	
	famFile_chr = system(paste0("ls gwasImputed_chr*.fam"), intern=TRUE)	
	bfile_chr = paste0(bedFile_chr, " ", bimFile_chr, " ", famFile_chr)
	fA = paste0(gsub(".bed", "", bedFile_chr[1]))
	tobeMerged = bfile_chr[-1]
	mergefilesetname = paste0("mergeGwasImputed.txt")
	write.table(tobeMerged, file=mergefilesetname, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	arg <- paste0(plink, " --bfile ", fA, " --merge-list ", mergefilesetname, " --make-bed --out gwasImputed")
	system(arg) 
	setwd('..')

} 

 
 







 

##########################################################################
########################################################################## filterImputeData.R
# Description:
# organize the files into the right folder and filter variants with relatively good quality. 

# Arguments:
# infoScore: the score of imputation quality for each variant.

library(doParallel) 
filterImputeData <- function(plink, infoScore){ 

	system("mv ./6-postImpute/gwasImputed* ./7-finalResults/") 

	 #####  get all snp info >>  
	# read each .impute2_info file, remove 1st line, add to another file and repeat  
	setwd("./5-imputeResults")
	## get all impute2_info files for each chunk
	files = system("ls *.impute2_info ", intern=TRUE)
	for (i in 1:length(files)){ 
	 	system(paste0("sed 1d ", files[i], "  >> impute2info.txt"))
	}  
	 
	### only SNPs and with two alleles 
	system("grep 'rs' impute2info.txt | awk '{if(length($4)==1 && length($5)==1) print}' | awk '{print $2, $7}' > impute2infoUpdate.txt") 
	## change dir 
	setwd("..") 
	## mv impute2info with added colnames and also filter impute2infoFiltered to 6-postImpute
	impute2info <- read.table(file="./5-imputeResults/impute2infoUpdate.txt", stringsAsFactors=F)  
	colnames(impute2info) = c("rs_id", "info") 
	 
	 ###############  filtering   
	snpWithBadInfo = impute2info[which(impute2info[, "info"] < infoScore), 1] 
	setwd("./7-finalResults")
	snpWithBadInfofn = "snpWithBadInfo.txt"
	write.table(snpWithBadInfo, file=snpWithBadInfofn, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=" ")
	### extract filtered SNPs  
	system( paste0(plink, " --bfile gwasImputed --exclude ", snpWithBadInfofn, " --make-bed --out gwasImputedFiltered")) 
	setwd("..")  # 
	setwd("..") ## /4-imputation

}


