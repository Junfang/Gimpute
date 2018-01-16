# File   : outputFinal.R
# Author : Junfang Chen
# Version: 15 AUG 2016

## Linking the output files
 
## imputed dataset
system( paste0("scp ./1-conversion/1_01_metaData.txt ./6-finalResults/metaData.txt " ) )
system( paste0("scp ./2-QC/2_01_snpHetXInstNumberAfter.txt ./6-finalResults/snpHetXInstNumberAfter.txt " ) )

system( paste0("scp ./2-QC/2_12_snpHweValuesAutoCt.txt ./6-finalResults/snpHweValuesAutoCt.txt " ) )
system( paste0("scp ./2-QC/2_13_snpHweValuesAlloCt.txt ./6-finalResults/snpHweValuesAlloCt.txt " ) )

## not for datasets only with controls
system( paste0("scp ./2-QC/2_13_snpHweValuesAlloAll.txt ./6-finalResults/snpHweValuesAlloAll.txt " ) ) 
system( paste0("scp ./2-QC/2_14_snpHweValuesAutoAll.txt ./6-finalResults/snpHweValuesAutoAll.txt " ) )  
 


## go to the directory/ dataset name 
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.bed ./6-finalResults/imputedSnpsDataset.bed") )
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.bim ./6-finalResults/imputedSnpsDataset.bim") )
system( paste0("scp ./4-imputation/4_6_removedSnpMissPostImp.fam ./6-finalResults/imputedSnpsDataset.fam") )

system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.bed ./6-finalResults/specificSnpsDataset.bed") )
system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.bim ./6-finalResults/specificSnpsDataset.bim") )
system( paste0("scp ./5-reductAndExpand/5_4_extSpecificDiffPos.fam ./6-finalResults/specificSnpsDataset.fam") )



 
 

