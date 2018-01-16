 
# File   : reductAndExpand.R
# Author : Junfang Chen
# Version: 15 AUG 2016


## subset imputed data into smaller set.
## define the directory used for calling r source codes 
  
##  


  
inputPrefix = "4_6_removedSnpMissPostImp"
inputOriginalQCed = "3_1_liftedDataset"

reducedToSpecificfn = "5_1_reducedToSpecific"
extSpecificDiffAllelefn = "5_2_extSpecificDiffAllele"
extSpecificMissPosfn = "5_3_extSpecificMissPos"
extSpecificDiffPosfn = "5_4_extSpecificDiffPos" 
 

## go to the directory/ dataset name
## imputed dataset
system( paste0("scp ./4-imputation/", inputPrefix, ".bed ./5-reductAndExpand/ " ) )
system( paste0("scp ./4-imputation/", inputPrefix, ".fam ./5-reductAndExpand/ " ) )
system( paste0("scp ./4-imputation/", inputPrefix, ".bim ./5-reductAndExpand/ " ) )

## original but QC-ed dataset
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".bed ./5-reductAndExpand/ " ) )
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".fam ./5-reductAndExpand/ " ) )
system( paste0("scp ./3-lifting/", inputOriginalQCed, ".bim ./5-reductAndExpand/ " ) )
## 
system( paste0("scp ./3-lifting/3_4_snpImpRefAlleles.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_4_snpDiffAlleles.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_3_snpMissPos.txt ./5-reductAndExpand/ ") )
system( paste0("scp ./3-lifting/3_2_snpDiffNamePos.txt ./5-reductAndExpand/ ") )


setwd("5-reductAndExpand/")
# 1. Reduce the imputed dataset to the SNPs before imputation.
plinkB = "/data/noether/tools/plink/plink --bfile " 
system(paste0(plinkB, inputPrefix, " --extract 3_4_snpImpRefAlleles.txt --make-bed --out ", reducedToSpecificfn) ) 

# 2. Add the SNPs with different alleles with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_4_snpDiffAlleles.txt"))==0 ){  
	system(paste0("scp ", reducedToSpecificfn, ".bed ", extSpecificDiffAllelefn, ".bed") ) 
	system(paste0("scp ", reducedToSpecificfn, ".bim ", extSpecificDiffAllelefn, ".bim") ) 
	system(paste0("scp ", reducedToSpecificfn, ".fam ", extSpecificDiffAllelefn, ".fam") ) 
} else { 
	system(paste0(plinkB, inputOriginalQCed, " --extract 3_4_snpDiffAlleles.txt --make-bed --out tmp") ) 
	system(paste0(plinkB, reducedToSpecificfn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificDiffAllelefn) ) 
	system("rm tmp.*")
  } 
  	
# 3. Add the SNPs with missing positions with their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_3_snpMissPos.txt"))==0 ){  
	system(paste0("scp ", extSpecificDiffAllelefn, ".bed ", extSpecificMissPosfn, ".bed") ) 
	system(paste0("scp ", extSpecificDiffAllelefn, ".bim ", extSpecificMissPosfn, ".bim") ) 
	system(paste0("scp ", extSpecificDiffAllelefn, ".fam ", extSpecificMissPosfn, ".fam") ) 
} else {
	system(paste0(plinkB, inputOriginalQCed, " --extract 3_3_snpMissPos.txt --make-bed --out tmp") ) 
	system(paste0(plinkB, extSpecificDiffAllelefn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificMissPosfn) ) 
	system("rm tmp.*")
  }
  	
# 4. Add the SNPs with different positions by their values from the dataset before removing SNPs. 
if ( file.size(paste0("3_2_snpDiffNamePos.txt"))==0 ){  
	system(paste0("scp ", extSpecificMissPosfn, ".bed ", extSpecificDiffPosfn, ".bed") ) 
	system(paste0("scp ", extSpecificMissPosfn, ".bim ", extSpecificDiffPosfn, ".bim") ) 
	system(paste0("scp ", extSpecificMissPosfn, ".fam ", extSpecificDiffPosfn, ".fam") ) 
} else {
	system(paste0(plinkB, inputOriginalQCed, " --extract 3_2_snpDiffNamePos.txt --make-bed --out tmp") ) 
	system(paste0(plinkB, extSpecificMissPosfn, " --bmerge  tmp.bed tmp.bim tmp.fam --make-bed --out ", extSpecificDiffPosfn) ) 
	system("rm tmp.*")
  }
 

system( paste0("rm ", inputPrefix, "*") ) 
system( paste0("rm ", inputOriginalQCed, "*") ) 
 
 )
system( paste0("rm  *.txt *.log *.hh") ) 
setwd("..")
 

 