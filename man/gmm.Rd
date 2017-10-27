\name{MLM Trials}
\alias{gmm}
\title{
Genomic mixed model
}
\description{
This function was developed to solve a mixel model for multi-environmental trials and/or replicated trials when genomic is available. The model includes a semi-parametric term to account for spatial variation when the field layout of the experiment is know. Covariates (fixed effects), genetics (random effect) and spatial term are fitted all in a single step.
}
\usage{
gmm(y,gen,dta=NULL,it=75,bi=25,th=1,model="BRR",...)
}
\arguments{
  \item{y}{
Numeric vector of phenotypes of length \eqn{obs} (missing values are allowed). \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix, with dimension \eqn{n} x \eqn{p}. Attention: Rows must be named with genotype IDs. Missing values are replaced by the SNP expectation.
}
  \item{dta}{
Data frame with \eqn{obs} number of rows containg the genotype ID, spatial information and any other set covariates. Make sure to add a column called "ID" (with capital letters) informing the genotype ID with some match in the object \code{gen}. For the spatial adjustment, it is necessary to add three numeric columns (not factors) to \code{dta}: "Block", "Row" and "Col" (case sensitive), where Block may refer to a field block or an enviroment without blocks Therefore, make sure that blocks from different environments are named differently. Row and Col provide the coordinates for the identification of neighbor plots.
}
  \item{it}{
Integer. Total numeric of MCMC iterations used to fit the model.
}
  \item{bi}{
Integer. Burn-in of MCMC iterations, i.e., number of iteration to be discarted prior to model convergence.
}
  \item{th}{
Integer. Thinning parameter: saves only 1 interation every \code{th}. Thinning is used to reduce the auto-correlation between Markov chains.
}
  \item{model}{
Prediction model: The options are: \code{BRR}, \code{BayesA}, \code{GBLUP}, \code{RKHS} and \code{RF}.
}
  \item{...}{
Pass arguments to the function that builds the spatial splines \code{NNsrc}: rho and dist. By default, \code{rho=1} and \code{dist=3}. To check how it looks like in the field, type \code{NNsrc(rho=1,dist=3)}. If \code{model="RF"}, it also passes the setting of random forest: number of trees \code{(NTR=50)}, sample size \code{(SSZ=80)}, max nodes \code{(MNO=10)}, node size \code{(NSZ=4)}, number of variables \code{(MTR=50)}, importance sampling \code{(IMP=TRUE)} and replacement \code{(RPL=TRUE)}.
}

}
\details{

The general model is \eqn{y=Xb+Zu+f(x)+e}, where \eqn{y} is the response variable, \eqn{Xb} refers to the fixed effects, \eqn{Zu} regards the genetic effect, \eqn{f(x)} represents the field variation, and \eqn{e} is the vector of residuals. In this model \eqn{u} is a link function that represents the genetic component, which depends on the model specified.

For whole-genome regression models (BRR or BayesA), \eqn{u = Ma}, where \eqn{M} is the matrix of genotypes. For kernel models (RKHS and GBLUP), \eqn{u=N(0,K\sigma2a)}, where K is either a Gaussian kernel (RKHS) or a linear kernel (GBLUP). For the non-parametric model (RF), \eqn{u} is updated every round of MCMC and prediction on unobseved genotyped are performed in every round of MCMC and average out in the end. To avoid over-representation of genotypes, \eqn{u} is not weighted according to the number of observations of each genotype.

Unobserved genotypes not provided in \code{dta} but provided in \code{gen} are predicted in the output of the function. Genotypes without genotypic information are transfered to the fixed effect (eg. checks). Missing loci are imputed with the expectation. If \code{dta} is not provided, the function will work as a regular genomic prediction model, so the length of \code{y} must match the number of rows of \code{gen}.

In whole-genome regression models, the regularization of the genetic term is either based on chosen prior (t, Gaussian), Gaussian (from ridge regression) and t (from BayesA). Kernel models (GBLUP and RKHS) are regularized as Gaussian process, which is similar to a ridge regression of Eigenvectors where the regularization of Eigenpairs also relies on the Eigenvalues. Random forest is, at some extend, regularized by multiple parameters including sample size (SSZ), numer of trees (NTR), number of variables (MTR), etc.

If there is a large number of trials and users acknowledge the necessity of sparse matrices, we recommend installing the Matrix package and run the following code that enables sparsity:

\code{source(system.file("add","sparseGMM.R",package="NAM"))}

}
\value{
The function gmm returns a list containing the fitted values (\code{hat}), observed values (\code{obs}), intercept (\code{mu}, incidence matrix of genotypes (\code{Z}) and the vector of breeding values (\code{EBV}). If fixed effects are provided, it also returns the design matrices and coefficients of fixed effects (\code{X},\code{b}). If the model was kernel or regression, output will include
the random effect coefficients (\code{g}), variance components of markers (\code{Vg}) and residuals (\code{Ve}). Kernel models regress the Eigenvectors of the kernel, weighted by the Eigenvalue. The coefficient (\code{cxx}) used in the \code{BRR} model to convert marker variance \eqn{Vb} into genetic variance \eqn{Va}. If spatial information is provided, the output includes the fitted spatial term (\code{sp}).
}

\author{
Alencar Xavier
}
\examples{

# For a demo with wulti-environmental trial type:
# demo(fittingMET)

# Checking heritability
data(tpod)
fit = gmm(y,gen,model = 'BRR')
fit$Vg*fit$cxx / (fit$Vg*fit$cxx+fit$Ve)

}
\keyword{ BLUP }
\keyword{ BGS }