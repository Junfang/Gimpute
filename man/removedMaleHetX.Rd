% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeQC.R
\name{removedMaleHetX}
\alias{removedMaleHetX}
\title{Remove male subjects with haploid heterozygous SNPs}
\usage{
removedMaleHetX(plink, inputPrefix, hhSubjCutOff, outputPrefix,
  outputFile_subjHetFreqAll, outputFile_subjHetFreqRetained,
  outputFile_SNPhhFreqAll)
}
\arguments{
\item{plink}{an executable PLINK program in either the current working directory or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK format files.}

\item{hhSubjCutOff}{the cutoff for removing male subjects with haploid heterozygous SNPs on the chromosome X.}

\item{outputPrefix}{the prefix of the output PLINK format files.}

\item{outputFile_subjHetFreqAll}{the output pure text file that stores male subjects that have heterozygous SNPs with their frequency (if any), i.e. the number of .hh SNPs in this male.}

\item{outputFile_subjHetFreqRetained}{the output pure text file that stores male subjects that have heterozygous SNPs with their frequency after 'improper' subject removal
(if any).}

\item{outputFile_SNPhhFreqAll}{the output pure text file that stores all heterozygous SNPs with their frequency (the number of males for this SNP), if any.}
}
\value{
The output PLINK format files and two pure text files (if any) with heterozygous SNPs and their respective frequency.
}
\description{
Determine the frequency of
}
