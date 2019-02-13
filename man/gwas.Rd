\name{GWAS}
\alias{gwas}
\alias{gwas2}
\alias{gwas3}
\alias{gwasGE}
\alias{meta3}
\title{
Empirical Bayes Genome Wide Association Mapping 
}
\description{
The gwas function calculates the likelihood ratio for each marker under the empirical Bayesian framework. The method allows analysis with multiple populations. \code{gwas2} is computationally optimized. \code{gwas3} was design for multiple random populations.
}
\usage{
gwas(y,gen,fam=NULL,chr=NULL,window=NULL,fixed=FALSE)
gwas2(y,gen,fam=NULL,chr=NULL,fixed=FALSE,EIG=NULL,cov=NULL)
gwas3(y,gen,fam=NULL,chr=NULL,EIG=NULL,cov=NULL)
gwasGE(Phe,gen,fam,chr=NULL,cov=NULL,ge=FALSE,ammi=1)
meta3(ByEnv,ammi=1)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozygous, heterozygous and reference homozigous. \code{NA} is allowed.
}
  \item{fam}{
Numeric vector of length \eqn{n} indicating a stratification factor or which subpopulation (\eqn{e.g.} family) that each observation comes from. Default assumes that all observations are from the same populations.
}
  \item{chr}{
Numeric vector indicating the number of markers in each chromosome. The sum of \eqn{chr} must be equal to the number of columns in \eqn{gen}. Default assumes that all markers are from the same chromosome.
}
  \item{window}{
Numeric. If specified, genetic distance between markers is used for moving window strategy. Window must be specified in Morgans (\eqn{e.g.} 0.05 would represent 5cM). Genetic distance is calculated assuming that individuals are RILs.
}
  \item{fixed}{
Logical. If TRUE, markers are treated as fixed effect and hence, evaluated through Wald statistics. If markers are specief as fixed, the argument 'window' is not applicable.
}
  \item{EIG}{
Output of the R function 'eigen'. It is used for user-defined kinship matrix.
}
  \item{cov}{
Numeric vector of length \eqn{n} to be used as covariate in the association analysis.
}
  \item{Phe}{
Numeric matrix of observations (\eqn{n*e}) where rows represent genotypes and columns represent environments. \code{NA} is allowed.
}
  \item{ge}{
Logical. If TRUE, meta-analysis (function gwasGE) will be done for the \eqn{G}x\eqn{E} interactions term only. If FALSE, variance components will be computed for three terms: genotype, environment and interaction.
}
  \item{ammi}{
Integer. It indicates the number of principal components used to represent \eqn{G}x\eqn{E} interactions through additive main-effects and multiplicative interaction (AMMI).
}
  \item{ByEnv}{
List of objected output from \code{gwas3} to perform meta-analysis.
}
}
\details{

Empirical Bayes model (Wang 2016) with a special incidence matrix is recreated to optimize the information provided by the subpopulations. Each locus is recoded as a vector with length \eqn{f} equal to number of subpopulations, or NAM families, as the interaction locus by family. For example, a locus heterozigous from an individual from subpopulation 2 is coded as [ 1, 0, 1 , ... ,\eqn{f} ], a locus homozigous for the reference allele from any subpopulation is coded as [ 2, 0, 0, ... , \eqn{f} ] and a locus homozigous for the founder allele from an individual from subpopulation 1 is coded as [ 0, 2, 0, ... ,\eqn{f} ]. The base model for genome scanning is described by:

\deqn{y = Xb + Zu + g + e}

That includes the fixed effect (\eqn{Xb}), the marker (\eqn{Zu}), the polygene (\eqn{g}) and the residuals (\eqn{e}). If the \eqn{window} term is specified, the model for genome scanning is expanded as follows:

\deqn{y = Xb + Zu[k-1] + Zu + Zu[k+1] + g - g[k] + e}

This model includes three extra terms: the left side genome ( \eqn{Zu[k-1]} ) and the right side genome ( \eqn{Zu[k+1]} ), also subtracting the window polygene ( \eqn{-g[k]} ). Windows are based on genetic distance, which is computed using Kosambi map function. The recombination rate is estimated under the assuption markers are ordered and that genotypes are recombinant inbred lines.

The polygenic term is calculated only once (Zhang et al 2010) using eigendecomposition with a GEMMA-like algorithm (Zhou ans Stephens 2012). Efficient inversion of capacitance matrix is obtained through the Woodbury matrix identities. Models and algorithms are described with more detail by Xavier et al (2015) and Wei and Xu (2016).

In order to analyze large dateset, one can avoid memory issues by using the function \eqn{gwas2}, but that the argument 'window' is not implemented for \eqn{gwas2}. This function also allows used-defined kindship through the argument EIG, and the use of a numeric covariate vector through the argument cov.

When multi-environmental trials are the target of mapping, one may use the function \eqn{gwasGE} to perform analysis by environment, followed by "meta-analysis" used to combine the results. This strategy provides an idea of the variation on QTL effect due to environment, genetic background (provided by the stratification factor) and the interaction between environment and genetics.

An alternative to this method is the mega-analysis, where one can provide the stratification factor as a combination of subpopulation and environment. Meta-analysis can be performed in a single step with function \code{gwasGE}, or users can perform multiple association analyses using \code{gwas3} and perform meta-analysis with \code{meta3}. In \code{gwasGE}, the same genotype will often appear more than once in the phenotypic and genotypic data, so that phenotypes are provided as a matrix. The statistical detail about the meta-analysis are available in the vignette \eqn{Background for Meta-analysis}.

The function \code{gwas3} is an alternative for association analysis and meta-analysis, also solved in the Empirical-Bayes framework for multiple populations. Unlike \code{gwas}, \code{gwas2} and \code{gwasGE}, this function does not set a reference allele and analysis each marker as the interaction of allele by stratification factor (ie. family or subpopulation). Therefore,  \code{gwas3} is compatible with any allele coding.

For further statistical background:

1) \code{system(paste('open',system.file("doc","gwa_description.pdf",package="NAM")))}

2) \code{system(paste('open',system.file("doc","gwa_ge_interactions.pdf",package="NAM")))}

}
\value{
The function gwas returns a list containing the method deployed (\eqn{Method}), a summary of predicted parameters and statistical tests (\eqn{PolyTest}), estimated genetic map for NAM panels (\eqn{MAP}) and the marker names (\eqn{SNPs}).
}
\references{

Wang, Q., Wei, J., Pan, Y., & Xu, S. (2016). An efficient empirical Bayes method for genomewide association studies. Journal of Animal Breeding and Genetics, 133(4), 253-263.  
  
Wei, J., & Xu, S. (2016). A Random Model Approach to QTL Mapping in Multi-parent Advanced Generation Inter-cross (MAGIC) Populations. Genetics, 202(2), 471-486.

Xavier, A., Xu, S., Muir, W. M., & Rainey, K. M. (2015). NAM: Association Studies in Multiple Populations. Bioinformatics, 31(23), 3862-3864.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.

Zhou, X., & Stephens, M. (2012). Genome-wide efficient mixed-model analysis for association studies. Nature genetics, 44(7), 821-824.
}
\author{
Alencar Xavier, Tiago Pimenta, Qishan Wang and Shizhong Xu
}
\examples{
\dontrun{
data(tpod)
gen=reference(gen)
gwa=gwas2(y=y,gen=gen,fam=fam,chr=chr,fixed=TRUE)
plot(gwa,pch=20,lwd=4)
}
}
\keyword{ NAM }
\keyword{ GWAS }
\keyword{ BLUP }
