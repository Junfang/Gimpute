% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phaseImpute.R
\name{prePhasingByShapeit}
\alias{prePhasingByShapeit}
\title{Prephasing genotypes using SHAPEIT}
\usage{
prePhasingByShapeit(shapeit, chrs, dataDIR, prefix4plinkEachChr, impRefDIR,
  phaseDIR, nThread, effectiveSize, nCore)
}
\arguments{
\item{shapeit}{an executable SHAPEIT program in either the current working directory or somewhere in the command path.}

\item{chrs}{specifiy the chromosomes for phasing.}

\item{dataDIR}{the directory where genotype PLINK files are located.}

\item{prefix4plinkEachChr}{the prefix of PLINK files for each chromosome.}

\item{impRefDIR}{the directory where the imputation reference files are located.}

\item{phaseDIR}{the directory where resulting pre-phased files will be located.}

\item{nThread}{the number of threads used for computation.}

\item{effectiveSize}{this parameter controls the effective population size.}

\item{nCore}{the number of cores used for computation. This can be tuned along with nThread.}
}
\value{
The pre-phased haplotypes for given chromosomes.
}
\description{
Perform prephasing for study genotypes by SHAPEIT for the autosomal and sex chromosome haplotypes using a reference panel (pre-set).
If ChrX is available then it is done differently by passing the flag --chrX to SHAPEIT.
}
\author{
Junfang Chen <junfang.chen3@gmail.com>
}
