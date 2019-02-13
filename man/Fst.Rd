\name{FST}
\alias{Fst}
\alias{plot.fst}
\title{
  Fixation Index
}
\description{
Genetic variation associated with markers distributed among subpopulations. The function generates a plot for structure diagnosis.
}
\usage{
  Fst(gen,fam)
}
\arguments{
  \item{gen}{
    Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers. SNPs must be coded as \code{0, 1, 2}, for founder homozigous, heterozygous and reference homozygous. \code{NA} is allowed.
  }
  \item{fam}{
    Numeric vector of length (\eqn{n}) indicating which subpopulations (\eqn{i.e.} family) each observation comes from. \code{NA} is not allowed.
  }
}
\details{
F-statistics (Wright 1965) represent the differentiation among populations for a given locus. Weir and Cockerham (1984) provided an unbiased version for molecular analysis. 

FIT is the correlation between gametes that unite to produce the individuals, relative to the gametes of the total population. FIS is the average over all subdivisions of the correlation between uniting gametes relative to those of their own subdivision. FST is the correlation between random gametes within subdivisions, relative to gametes of the total population. Neutral markers have an expected FST 0.05.
}
\value{
List with values of FST, FIS and FIT. Unbiased F-statistics from weighted AOV (Weir and Cockerham 1984).
}
\references{

Weir, B. S., and Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. Evolution, 38(6), 1358-1370.

Wright, S. (1965). The interpretation of population structure by F-statistics with special regard to systems of mating. Evolution, 19(3), 395-420.

}
\author{
  Alencar Xavier and William Muir
}
\examples{
  \dontrun{
  data(tpod)
  Fstat = Fst(gen=gen,fam=fam)
  plot(Fstat,chr=chr)
    }
}
\keyword{ Fst }
