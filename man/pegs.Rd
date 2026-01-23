\name{MLM MV-PEGS}
\alias{pegs}
\alias{Y}
\alias{X}
\alias{simu}
\title{Pseudo-Expectation Gauss-Seidel for Multivariate Models}
\usage{
pegs(Y, X, maxit = 100, logtol = -4, NNC = FALSE, 
     covbend = 1.1, covMinEv = 10e-4, XFA = -1)
}
\arguments{
\item{Y}{A numeric matrix (n x k) of dependent variables, where n is the number of
observations and k is the number of response variables. Missing values should be NA.}

\item{X}{A numeric matrix (n x p) of predictors (i.e., markers) for a single random
effect, or a named list of such matrices for multiple random effects.}

\item{maxit}{An integer specifying the maximum number of iterations.}

\item{logtol}{A numeric value for the convergence tolerance on a log10 scale.
The algorithm stops when the log10 of the sum of squared changes in
coefficients is less than this value.}

\item{NNC}{A logical value. If TRUE, imposes a constraint of non-negative correlations
on the genetic covariance matrix.}

\item{covbend}{A numeric bending factor used to inflate the diagonal of the covariance
matrix if it is not positive semi-definite. The bending amount is calculated as
\code{abs(min_eigenvalue) * covbend}.}

\item{covMinEv}{The minimum allowable eigenvalue for the genetic covariance matrix.
If the smallest eigenvalue falls below this threshold, the matrix is "bent"
(diagonally inflated) to make it positive semi-definite.}

\item{XFA}{An integer specifying the number of factors (principal components) to use in
the factor-analytic approximation of the genetic covariance matrix.
\itemize{
  \item \code{XFA = -1} (Default): Fits a full-rank model (no approximation).
  \item \code{XFA = 0}: Fits a diagonal model, assuming independent genetic effects among traits.
  \item \code{XFA > 0}: Fits a reduced-rank model using the specified number of factors.
}}
}
\value{
A list containing the following components:
\describe{
  \item{mu}{Vector of intercepts for each response variable.}
  \item{b}{A matrix of regression coefficients (p x k) if a single \code{X} matrix is provided,
    or a named list of such matrices if a list is provided for \code{X}.}
  \item{hat}{Matrix of predicted values (n x k).}
  \item{h2}{Vector of heritability estimates for each trait.}
  \item{GC}{The genetic correlation matrix (k x k) if a single \code{X} matrix is provided,
    or a named list of such matrices for each random effect.}
  \item{bend}{The inflation factor added to the diagonal of the genetic variance matrix.
    This will be a single value if one random effect is fit, or a named vector if multiple
    effects are fit.}
  \item{numit}{The number of iterations until convergence.}
  \item{cnv}{The final convergence value.}
}
}
\description{
This function implements the Pseudo-Expectation Gauss-Seidel (PEGS) algorithm for
multi-trait genomic prediction models. It supports fitting one or more random effects and
can apply factor-analytic approximations to the genetic covariance structure.
}
\examples{
# Load simulated data with one random effect
data(test)

# Fit a standard multi-trait model
fit <- pegs(Y, X)

# Heritability
print(fit$h2)

# Genetic correlations
print(fit$GC)

# Accuracy
print(mean(diag(cor(fit$hat, simu$tbv))))

}
\references{

Xavier, A., Runcie, D., & Habier, D. (2025). Megavariate methods capture complex genotype-by-environment interactions. Genetics, 229(4), iyae179.

Xavier, A., & Habier, D. (2022). A new approach fits multivariate genomic prediction models efficiently. Genetics Selection Evolution, 54(1), 45.

}