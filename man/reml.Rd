\name{MLM REML}
\alias{reml}
\alias{MCreml}
\title{
Restricted Maximum Likelihood
}
\description{
Univariate REML estimators and variance components for a single random variable fitted by an EMMA-like algorithm.
}
\usage{
reml(y,X=NULL,Z=NULL,K=NULL)
MCreml(y,K,X=NULL,MC=300,samp=300)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Formula or incidence matrix (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{Z}{
Formula or numeric matrix (\eqn{n} by \eqn{p}) that corresponds to the incidence matrix of random effect. \code{NA} is not allowed.
}
  \item{K}{
Numeric matrix (\eqn{p} by \eqn{p}). Kinship matrix for random effect with \eqn{p} parameters. \code{NA} is not allowed.
}
  \item{MC}{
Number of sampling procedures to estimate variance components using MCreml.
}
  \item{samp}{
Sample size of the sampling procedure to estimate variance components using MCreml.
}
}
\details{
Solve mixed models with a single random effects minizing the log restricted maximum likelihood (REML) using the EMMA algorithm (Kang et al 2008). Prediction of random coefficients for ridge-type model are performed according to VanRaden (2008), and kernel-type model via RKHS according to de los Campos et al. (2010).

If y is a matrix with multiple traits, the fuctions solves the mixed model via an ECM algorithm adapted from the EMMREML package (Akdemir and Godfrey 2014).

MCreml is based on subsampling with \code{samp} observations at a time, repeating the procedure \code{MC} times. Subsampling is a common Monte Carlo strategy to reduce the computational burden to estimate variance components in large datasets.

}
\value{
The function reml returns a list with variance components and heritability (VC), fixed effect coefficients and standard variations (Fixed) and estimated breeding values (EBV). 
}
\references{

Akdemir, D., and O. U. Godfrey (2014) EMMREML: Fitting Mixed Models with Known Covariance Structures. R Package Version 2.0. Available at: http://CRAN.R-project.org/package=EMMREML.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D., Daly, M. J., & Eskin, E. (2008). Efficient control of population structure in model organism association mapping. Genetics, 178(3), 1709-1723.

VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), 4414-4423.

}
\author{
Alencar Xavier, Tiago Pimenta and Shizhong Xu
}
\examples{

# Fitting a random model
data(tpod)
FIT = reml(y=y,Z=~as.factor(fam))

# Fitting GBLUP
G = GRM(gen)
GBLUP = reml(y=y,K=G)

# GBLUP vs RRBLUP
g = tcrossprod(gen)
gblup = reml(y=y,K=g)
rrblup = reml(y=y,Z=gen)
rbind(gblup$VC,rrblup$VC)
gebv_gblup = gblup$EBV
gebv_rrblup = c(tcrossprod(t(rrblup$EBV),gen))
plot(gebv_gblup,gebv_rrblup)

}
\keyword{ REML }
\keyword{ BLUP }
