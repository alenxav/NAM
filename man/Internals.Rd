\name{Optimizations}
\alias{calcSize}
\alias{funI}
\alias{funX}
\alias{inputRow}
\alias{RcppExports}
\alias{timesMatrix}
\alias{timesVec}
\alias{SAMP}
\alias{SAMP2}
\alias{gs}
\alias{KMUP}
\alias{KMUP2}
\alias{emRR}
\alias{emBA}
\alias{emBB}
\alias{emBC}
\alias{emBL}
\alias{emDE}
\alias{emEN}
\alias{emML}
\alias{emML2}
\alias{emMX}
\alias{emCV}
\alias{BRR2}
\alias{Gdist}
\alias{markov}
\alias{eigX}
\alias{covar}
\alias{PedMat}
\alias{PedMat2}
\alias{LD}
\alias{G2A_Kernels}
\alias{Import_data}
\alias{NNsrc}
\alias{NNcov}
\alias{CNT}
\alias{MSX}
\alias{IMP}
\alias{NOR}
\alias{GAU}
\alias{GRM}
\alias{SPC}
\alias{SPM}
\alias{press}
\alias{emGWA}
\alias{BCpi}
\alias{mrr}
\title{
  Internal functions
}
\description{
Internal function under optimization, complimentary statistics, and loops written in C++ to speed up \eqn{gwas}, \eqn{gibbs} and \eqn{wgr}.
}
\details{

Some of the functions available for users include:

01) \code{Import_data(file,type=c('GBS','HapMap','VCF')):}{  This function can be used to import genotypic data in the NAM format, providing a list with a genotypic matrix \code{gen} coded as \code{012} and a vector \code{chr} with count of markers per chromosome. Currently, it helps users to import three types of files: GBS text, HapMap and VCF.}

02) \code{markov(gen,chr):}{  Imputation method based forwards Markov model for SNP data coded as \code{012}. We recommend users to remove non-segregating markers before using this function.}

03) \code{LD(gen):}{  Computes the linkage disequilibrium in terms of r2 for SNP data coded as \code{012}. Missing data is not allowed.}

04) \code{PedMat(ped):}{  Builds a kinship from a pedigree. Input format is provided with \code{PedMat()}.}

05) \code{PedMat2(ped,gen=NULL,IgnoreInbr=FALSE,PureLines=FALSE):}{  Builds a kinship from a genomic data and pedigree. Useful when not all individuals are genotyped. Row names of \code{gen} must indicate the genotype id.}

06) \code{Gdist(gen, method = 1):}{  Computes genetic distance among individuals. Five methods are available: 1) Nei distance; 2) Edwards distance; 3) Reynolds distance; 4) Rogers distance; 5) Provesti's distance. 6) Modified Rogers distance}

07) \code{covar(sp=NULL,rho=3.5,type=1,dist=2.5):}{  Builds a spatial kernel from field plot information. Input format is provided with \code{covar()}. Parameter \code{rho} detemines the decay of relationship among neighbor plots. \code{type} defines if the kernel is exponential (1), Gaussian (2) or some intermediate. \code{dist} informs the distance ratio between range neighbors and row neighbors.}

08) \code{eigX(gen,fam):}{  Computes the input of the argument \code{EIG} of the function \code{gwas2}.}

09) \code{G2A_Kernels(gen):}{  Computes a list of orthogonal kernels containing additive, dominant and first-order epistatic effects, in accordance to the G2A model from ZB Zeng et al. 2005. These kernels can be used for description of genetic architecture through variance components, for that we recommend packages \code{varComp} and \code{BGLR}.}

10) \code{NNsrc(sp=NULL,rho=1,dist=3):}{  Using the same field data input required by the function \code{covar}, this function provides a list of nearest neighbor plots for each entry.}

11) \code{NNcov(NN,y):}{  This function utilizes the output of \code{NNsrc} to generate a numeric vector, averageing the observed values of \code{y}. This function is useful to generate field covariates to control micro-environmental variance without krigging.}

11) \code{emXX(y,gen,...):}{  Fits whole-genome regressions using the expectation-maximization algorithm as opposed to MCMC. Currently avaible methods include BayesA (\code{emBA}), BayesB (\code{emBB}), BayesC (\code{emBC}), BLASSO (\code{emBL}), BLASSO2 (\code{emDE}), Elastic-Net (\code{emEN}), maximum likelihood (\code{emML}), Ridge regression (\code{emRR}), and mixed-normals prior (\code{emMX}). A cross-validation option is also available (\code{emCV}).}

12) \code{CNT(X):}{  Centralizes parameters from matrix \code{X}.}

13) \code{IMP(X):}{  Imputes missing points from matrix \code{X} with the average value of the column.}

14) \code{GAU(X):}{ Creates a Gaussian kernel from matrix \code{X}.}

15) \code{GRM(X, Code012=FALSE):}{ Creates genomic relationship matrix as linear kernel from matrix \code{X}. If genotypes are coded as 012 and Code012=TRUE, the kinship is the same as proposed by VanRaden (2008), otherwise the outcome is an additive G2A kernel.}

16) \code{MSX(X):}{ Computes the cross-product of each column of \code{X} and the sum of variances of each column of \code{X}.}

17) \code{NOR(y,X,cxx,xx,maxit=50,tol=10e-6):}{ Solves a ridge regression using GSRU, where  \code{y} corresponds to the response variable, \code{X} is the set of parameters, \code{cxx} and \code{xx} are the output from the MSX function, \code{maxit} and \code{tol} are the convergence criteria.}

18) \code{SPC(y,blk,row,col,rN=3,cN=1):}{ Computes a spatial covariate, similar to what could be obtained using NNsrc and NNcov but in a single step. It often is faster than NNsrc/NNcov.}

19) \code{SPM(blk,row,col,rN=3,cN=1):}{ Computes a spatial matrix that capture nearest neighbots, to be used as design matrix of random effects. The least-square solution gives the same as SPC.}

20) \code{BRR2(y,X1,X2,it=1500,bi=500,df=5,R2=0.5):}{ A simple C++ implementation of a Bayesian Ridge Regression that accomodates two random effects.}

21) \code{emML2(y,X1,X2,D1=NULL,D2=NULL):}{ A simple C++ implementation of emML that accomodates two random effects.}

22) \code{press(y,K,MaxIt=10):}{Solves a PRESS-regularized GBLUP. You can provide K as a matrix or as the output of the functin eigen. MaxIt the maximum number of iterations to for updating missing values (if any) if H*y does not converge.} 

23) \code{emGWA(y,gen):}{ A vanilla algorithm written in C++ for GWAS (very simple, but very efficient). It fits a snpBLUP via EM-REML based GSRU, then run an additional round checkinkg the likelihood of treating each marker as fixed effect instead of random, thus avoiding double-fitting. It returns the marker p-values, snpBLUP marker effects for genomic prediction, LS marker effects from the GWAS, variance components, heritability, and GEBVs (fitted values).}

24) \code{BCpi(y,X,it=3000,bi=500,df=5,R2=0.5):}{ A vanilla implementation in C++ of BayesCpi for GWAS or GWP. It returns the marker p-values (as the minus log probability of marker excluded), marker effects for genomic prediction, probability of marker included, variance components, heritability, and GEBVs (fitted values). }

25) \code{mrr(Y,X):}{A C++ implementation for multivariate regression.}

}
\author{
  Alencar Xavier
}
\examples{
\dontrun{


# Forward gen imputation
data(tpod)
fast.impute = markov(gen,chr)

# Wright's A matrix 
PedMat()

# Pairwise LD
ld = LD(gen[,1:10])
heatmap(ld)

# Spatial correlation (kernel-based)
covar()

# Spatial correlation (NN-based)
NNsrc()

# Genetic distance
round(Gdist(gen[1:10,],method=1),2)

# PCs of a NAM kinship
eG = eigX(gen,fam)
plot(eG[[2]],col=fam)

# Polygenic kinship matrices
Ks = G2A_Kernels(gen)
ls(Ks)

# Genomic regression fitted via EM
h = emBA(y,gen)
plot(h$b,pch=20)

# GBLUP and RRBLUP
g = GRM(gen)
eg = eigen(g)
gblup = emML(y=y, gen=eg$vectors,D=eg$values)
rrblup = emML(y=y, gen=gen)
plot(gblup$hat,rrblup$hat,xlab = 'gblup',ylab='rrblup')

# Vanilla GWAS
gwa = emGWA(y,gen)
plot(gwa$PVAL,pch=20)
}
}
\keyword{ NAM }
\keyword{ BLUP }
\keyword{ GWAS }
\keyword{ LD }
\keyword{ GRM }
\keyword{ SPATIAL }
\keyword{ PEDIGREE }
