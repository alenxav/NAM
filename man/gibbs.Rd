\name{MLM Gibbs}
\alias{gibbs}
\alias{ml}
\alias{gibbs2}
\alias{plot.gibbs}
\title{
Bayesian Mixed Model
}
\description{
Mixed model solver through Bayesian Gibbs Sampling or iterative solution.
}
\usage{
gibbs(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,Iter=1500,Burn=500,
      Thin=2,DF=5,S=0.5,nor=TRUE,GSRU=FALSE)
ml(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,DF=5,S=0.5,nor=TRUE)
gibbs2(Y,Z=NULL,X=NULL,iK=NULL,Iter=150,Burn=50,Thin=1,DF=5,S=0.5,nor=TRUE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{Z}{
Right-hand side formula or list of numeric matrices (\eqn{n} by \eqn{p}) with incidence matrices for random effect. \code{NA} is not allowed.
}
  \item{X}{
Right-hand side formula or incidence matrices (\eqn{n} by \eqn{p}) for fixed effect. \code{NA} is not allowed.
}
  \item{iK}{
Numeric matrix or list of numeric matrices (\eqn{p} by \eqn{p}) corresponding to the the inverse kinship matrix of each random effect with \eqn{p} parameters.
}
  \item{iR}{
Numeric matrix (\eqn{n} by \eqn{n}) corresponding to the inverse residual correlation matrix.
}
  \item{Iter}{
Integer. Number of iterations or samples to be generated.
}
  \item{Burn}{
Integer. Number of iterations or samples to be discarted.
}
  \item{Thin}{
Integer. Thinning parameter, used to save memory by storing only one every 'Thin' samples.
}
  \item{DF}{
Integer. Hyperprior degrees of freedom of variance components.
}
  \item{S}{
Integer or NULL. Hyperprior shape of variance components. If NULL, the hyperprior solution for the scale parameter is calculated as proposed by de los Campos et al. (2013).
}
  \item{nor}{
Logical. If TRUE, it normilizes the response variable(s).
}
  \item{GSRU}{
Logical. If TRUE, it updates the regression coefficient using Gauss-Seidel Residual Update (Legarra and Misztal 2008). Useful for p>>n, but does not work when iK or iR are provided.
}
  \item{Y}{
Numeric matrix of observations for multivariate mixed models. Each column represents a trait. \code{NA} is allowed.
}
}
\details{
The general model is \eqn{y=Xb+Zu+e}, where \eqn{u=N(0,A\sigma2a)} and \eqn{e=N(0,R\sigma2e)}. The function solves Gaussian mixed models in the Bayesian framework as described by Garcia-Cortes and Sorensen (1996) and Sorensen and Gianola (2002) with conjugated priors. The alternative function, "ml", finds the solution iteratively using the full-conditional expectation. The function "gibbs2" can be used for the multivariate case, check Xavier et al. (2017) for an example of multivariate mixed model using Gibbs sampling.
}
\value{
The function gibbs returns a list with variance components distribution a posteriori (Posterior.VC) and mode estimated (VC.estimate), a list with the posterior distribution of regression coefficients (Posterior.Coef) and the posterior mean (Coef.estimate), and the fitted values using the mean (Fit.mean) of posterior coefficients.
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

Garcia-Cortes, L. A., and Sorensen, D. (1996). On a multivariate implementation of the Gibbs sampler. Genetics Selection Evolution, 28(1), 121-126.

Sorensen, D., & Gianola, D. (2002). Likelihood, Bayesian, and MCMC methods in quantitative genetics. Springer Science & Business Media.

Xavier, A., Hall, B., Casteel, S., Muir, W. and Rainey, K.M. (2017). Using unsupervised learning techniques to assess interactions among complex traits in soybeans. Euphytica, 213(8), p.200.

}
\author{
Alencar Xavier
}
\examples{

data(tpod)

# Fitting GBLUP
K = tcrossprod(gen-1)
K = K/mean(diag(K))
iK = chol2inv(K)
# FIT
test1 = gibbs(y,iK=iK,S=1)
# PLOT
par(mfrow=c(1,3))
plot(test1$Fit.mean,y,pch=20,lwd=2,col=3,main='GBLUP')
plot(test1,col=4,lwd=2)
# Heritability
print(paste('h2 =',round(test1$VC.estimate[1]/sum(test1$VC.estimate),3)))

# Fitting RKHS
E = as.matrix(dist(gen))
G = exp(-E^2/mean(E^2))
EIG = eigen(G,symmetric = TRUE)
ev = 20
U = EIG$vectors[,1:ev]
iV = diag(1/EIG$values[1:ev])
# FIT
test2 = gibbs(y,Z=U,iK=iV,S=1)
# PLOT
par(mfrow=c(1,3))
plot(test2$Fit.mean,y,pch=20,lwd=2,col=2,main='RKHS')
plot(test2,col=3,lwd=2)
# Heritability
print(paste('h2 =',round(test2$VC.estimate[1]/sum(test2$VC.estimate),3)))


}
\keyword{ BLUP }
\keyword{ BGS }