\name{SNP QC}
\alias{snpQC}
\alias{cleanREP}
\alias{reference}
\title{
SNP Quality Control
}
\description{
Functions for quality control. 'snpQC' may be used to count/remove neighbor repeated SNPs, markers with MAF lower than a given threshold, and imputations. 'cleanREP' identifies and merge duplicate genotypes. The 'reference' function changes the reference genotype. For NAM populations, this function must be used when genotypes are coded according to the reference genome instead of the standard parent.
}
\usage{
snpQC(gen,psy=1,MAF=0.05,misThr=0.8,remove=TRUE,impute=FALSE)
cleanREP(y,gen,fam=NULL,thr=0.95)
reference(gen,ref=NULL)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozigous, heterozigous and reference homozigous. \code{NA} is allowed.
}
  \item{psy}{
Tolerance parameter for markers in Perfect SYymmetry (psy). This QC remove identical markers (aka. full LD) that carry the same information. Default is 1, which removes only SNPs 100\% equal to its following neighbor.
}
  \item{MAF}{
Minor Allele Frequency. Default is 0.05. Useful to inform or remove markers below the MAF threshold. Markers with standard deviation below the MAF threshold will be also removed.
}
  \item{misThr}{
Missing value threshold. Default is 0.8, removing markers with more than 80 percent missing values. 
}
  \item{remove}{
Logical. Remove SNPs due to PSY or MAF.
}
  \item{impute}{
If TRUE, impute missing values using the expected value.
}
  \item{y}{
Numeric vector (\eqn{n}) or numeric matrix (\eqn{n} x \eqn{t}) of observations describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{fam}{
Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. Default assumes that all observations are from the same populations.
}
  \item{thr}{
Threshold above which genotypes are considered identical. Default is 0.95, merging genotypes >95 percent identical.
}
\item{ref}{
Numeric vector of length \eqn{n} with elements coded as \code{0, 1, 2}, it represents the genotypic information of a new reference genotype. Default assumes that more frequent allele represents the reference genome.
}
}
\value{
snpQC - Returns the genomic matrix without missing values, redundancy or low MAF markers.

cleanREP - List containing the inputs without replicates. Groups of replicates are replaced by a single observation with the phenotypic expected value. The algorithm keeps the genotypic information of the first individual (genotypic matrix order).

reference - Returns a recoded \eqn{gen} matrix
}
\author{
Alencar Xavier, Katy Rainey and William Muir
}
\examples{
  \dontrun{
data(tpod)
gen=reference(gen)
gen=snpQC(gen=gen,psy=1,MAF=0.05,remove=TRUE,impute=FALSE)
test=cleanREP(y,gen)
  }
}
