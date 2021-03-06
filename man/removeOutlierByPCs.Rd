% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeQC.R
\name{removeOutlierByPCs}
\alias{removeOutlierByPCs}
\title{Remove population outliers by principle components}
\usage{
removeOutlierByPCs(plink, gcta, inputPrefix, cutoff, cutoffSign,
  inputPC4subjFile, outputPC4outlierFile, outputPCplotFile, outputPrefix)
}
\arguments{
\item{plink}{an executable PLINK program in either the current working directory or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK format files.}

\item{cutoff}{the cutoff that would distinguish the outliers from ordinary population. 
If it is null, then no outliers or not necessary to remove outliers.
If it has only a value, then one also has to define the cutoffSign as below.
If it has two values, then one doesn't have to define the cutoffSign.}

\item{cutoffSign}{the cutoff sign: 'greater' or 'smaller' that would determine if the outliers should be greater or smaller than the cutoff value.
if the cutoff score has two values, then no need to define the cutoffSign.}

\item{inputPC4subjFile}{the pure text file that stores all the subject IDs after QC and their corresponding eigenvalues of the first two principle components, if any.}

\item{outputPC4outlierFile}{the pure text file that stores the outlier IDs and their corresponding eigenvalues of the first two principle components, if any.}

\item{outputPCplotFile}{the plot file for visualizing the first two principle components of all subjects without population outliers, if any.}

\item{outputPrefix}{the prefix of the output PLINK format files..}

\item{outputSNPs}{the output pure text file that stores the removed monomorphic SNPs.}
}
\value{
The output PLINK format files after outlier removal. The output pure text file (if any) for storing removed outlier IDs and their corresponding PCs. 
The plot file (if any) for visualizing the first two principle components after outlier removal.
}
\description{
This function decides if one would remove population outliers or not. If the outliers are necessary to be removed 
then one uses the eigenvalues from the first principle component as a criterion for finding out the outliers by defining a proper cutoff.
}
\author{
Junfang Chen <junfang.chen3@gmail.com>
}
