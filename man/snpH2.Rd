\name{SNP H2}
\alias{snpH2}
\alias{plot.H2}
\title{
SNP heritability
}
\description{
Calculates the ability of markers to carry a gene (Foneris et al. 2015). The index is also an indicator of Mendelian segregation.
}
\usage{
snpH2(gen,K=NULL)
}
\arguments{
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n} rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{K}{
Optional. Numeric matrix containing the genetic relationship matrix. A square and symmetric matrix (\eqn{n}x\eqn{n}) with the same observations as the matrix gen, also in the same order. If not provided, the kinship matrix is estimated from the genotypic matrix.
}
}
\value{
Numeric vector containing the heritability of each markers. Foneris et al. (2015) recommends to avoid using markers with index lower than 0.98. 
}
\references{
Forneris, N. S., Legarra, A., Vitezica, Z. G., Tsuruta, S., Aguilar, I., Misztal, I., & Cantet, R. J. (2015). Quality Control of Genotypes Using Heritability Estimates of Gene Content at the Marker. Genetics 199(3):675-681.
}
\author{

Alencar Xavier

}
\examples{

data(tpod)
Heritability=snpH2(gen)
plot(Heritability,chr=chr)

}