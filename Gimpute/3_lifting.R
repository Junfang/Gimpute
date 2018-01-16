
  
# File   : 3_lifting.R
# Author : Junfang Chen
# Version0: 28 Jun 2016
# VersionX: 15 Jan 2018
 	
  

 
## step 1  
system( paste0("cp ./2-QC/2_14_removedOutliersQCed.bed ./3-lifting/3_1_liftedDataset.bed") )
system( paste0("cp ./2-QC/2_14_removedOutliersQCed.fam ./3-lifting/3_1_liftedDataset.fam") )
system( paste0("cp ./2-QC/2_14_removedOutliersQCed.bim ./3-lifting/3_1_liftedDataset.bim") )
setwd("./3-lifting/")

   
## required functions
source(paste0(codeDir, "3_function4lifting.R"))

  

system( paste0("rm  *.log *.hh") ) 

setwd("..")

