% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snpInfoUpdate.R
\name{updatedSnpInfo}
\alias{updatedSnpInfo}
\title{Update the probes or SNPs information}
\usage{
updatedSnpInfo(plink, inputPrefix, chipAnnoFile, chipType, outputPrefix)
}
\arguments{
\item{plink}{an executable PLINK program in either the current working directory or somewhere in the command path.}

\item{inputPrefix}{the prefix of the input PLINK format files.}

\item{chipAnnoFile}{a pure text file that stores the chip annotation information (Affymetrix, Illumination, PsychChip and so on). 
This file can be found http://www.well.ox.ac.uk/~wrayner/strand/.}

\item{chipType}{a string name that defines the type of the chip annotation file: 'Illumina', 'affymetrix' or 'psychChip'.}

\item{outputPrefix}{the prefix of the output PLINK format files.}
}
\value{
The output PLINK format files after updating SNP information.
}
\description{
Remove SNP information including SNP name, genomic position, chromosome location and the strand information, with the help of the chip annotation file. 
This chip annotation file is defined in the configuration folder.
}
\author{
Junfang Chen <junfang.chen3@gmail.com>
}
