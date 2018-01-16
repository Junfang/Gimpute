 # 

# File   : align2ref.R
# Author : Junfang Chen
# Version0: 06 Aug 2016
# VersionX: 15 Jan 2018
 	
 

  
# Note it may be needed: dos2unix dataSetImpRefCompare.pl  
alignScript4compareName = paste0(codeDir,"utilities/dataSetImpRefCompareName.pl ")
alignScript = paste0(codeDir,"utilities/dataSetImpRefComparePos.pl ")
## dir
# impRefdir = "/data/noether/dataProcessResults/10_Common/imputeReference/ALL_1000G_phase1integrated_v3_impute_macGT1 "
impRefdir = "/data/noether/datatmp-nobackup/tb2refDatabase/imputeRef/1000Gphase1/ "
 
# 2. Remove SNPs for which the name (rs-ID) has a different position (i.e. combination of
# base pair position and chromosome) in the imputation reference files. 

inputPrefix = "3_1_liftedDataset" 
out2.snp = "3_2_snpDiffNamePos"
out2 = "3_2_removedSnpDiffNamePos"
 
inputBIMfn = paste0(inputPrefix, ".bim ") ## bim file from step1 
system( paste0(alignScript4compareName, inputBIMfn, impRefdir, " pos > ", out2.snp, ".txt"))
system( paste0(plink, " --bfile ",  inputPrefix, " --exclude ", out2.snp, ".txt --make-bed --out ", out2) )  


# 3. Remove SNPs which bp and chr position are not contained in the imputation reference files.

out3 = "3_3_removedSnpMissPos"
out3.snp = "3_3_snpMissPos"

 
bim2fn = paste0(out2, ".bim ")  
system( paste0(alignScript, bim2fn, impRefdir, " pos > ", out3.snp, ".txt"))
system( paste0(plink, " --bfile ",  out2, " --exclude ", paste0(out3.snp, ".txt"), " --make-bed --out ", out3) ) 



## 4. Remove SNPs which have an allele which is not in the imputation reference files for that SNP.
out4 = "3_4_removedSnpDiffAlleles"
out4.snp = "3_4_snpDiffAlleles"
out4.snpRetained = "3_4_snpImpRefAlleles"

 
bim3fn = paste0(out3, ".bim ")  
system( paste0(alignScript, bim3fn, impRefdir, " alleles > ", out4.snp, ".txt"))
system( paste0(plink, " --bfile ",  out3, " --exclude ", paste0(out4.snp, ".txt"), " --make-bed --out ", out4) )  

## 3.2 Remove SNPs which have an allele which is not in the imputation reference files for that SNP.
snpDifAllele = read.table(paste0(out4.snp, ".txt"), stringsAsFactors=F)
snpDifAllele = snpDifAllele[,1] 
inputBIM = read.table(paste0(out3, ".bim"), stringsAsFactors=F)
snpRetained = setdiff(inputBIM[,2], snpDifAllele)
write.table(snpRetained, file=paste0(out4.snpRetained, ".txt"), quote=F, row.names=F, col.names=F, eol="\r\n", sep=" ")


 
