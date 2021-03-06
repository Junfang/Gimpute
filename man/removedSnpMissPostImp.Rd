% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phaseImpute.R
\name{removedSnpMissPostImp}
\alias{removedSnpMissPostImp}
\title{Remove SNPs after post imputation}
\usage{
removedSnpMissPostImp(plink, inputPrefix, missCutoff, snpWithManyMissSNPfile,
  outputPrefix)
}
\arguments{
\item{plink}{an executable PLINK program in either the current working directory or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK format files.}

\item{missCutoff}{the cutoff of the least number of instances for a SNP that is not missing.}

\item{snpWithManyMissSNPfile}{the output file of SNPs with pre-defined missing values.}

\item{outputPrefix}{the prefix of the PLINK format files.}
}
\value{
The PLINK format files after post imputation quality control and a pure text file contains SNPs with pre-defined missing values.
}
\description{
Remove SNPs which have a non missing value for less than a predefined number of instances.
}
\examples{
  
}
\author{
Junfang Chen <junfang.chen3@gmail.com>
}
