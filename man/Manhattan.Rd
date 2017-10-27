\name{Manhattan}
\alias{plot.NAM}
\title{
Manhattan plot for Association Studies 
}
\description{
Generates a graphical visualization for the output of the function gwas/gwas2/gwasGE. 
}
\usage{
\method{plot}{NAM}(x,...,alpha=0.05,colA=2,colB=4,find=NULL,FDR=NULL,gtz=FALSE,phys=NULL)
}
\arguments{
   \item{x}{
Output of the gwas/gwas2 function.
}
  \item{...}{
Further arguments passed to or from other methods.
}
  \item{alpha}{
Numberic. Significance threshold to display in the Manhattan plot.
}
  \item{colA}{
Color of odd chromosomes in the Manhattan plot.
}
  \item{colB}{
Color of even chromosomes in the Manhattan plot.
}
  \item{find}{
Integer. If provided, you can click on the specified number of hits in the Manhattan plot to obtain the name of the markers.
}
  \item{FDR}{
Null or numeric between zero and one. If provided, it will display the Manhattan plot with Bonferroni threshold by chromosome, adjusted for the specified false discovery rate (FDR). Thus, zero provides the Bonferroni correction.
}
  \item{gtz}{
Logical. If TRUE, the argument FDR will just take into account markers with p-value Greater Than Zero (GTZ).
}
  \item{phys}{
Numeric vector with length equal to the number of markers. If provided, the Manhattan plot is generated using the physical distance in the x axis.
}
}
\author{
Alencar Xavier and William Beavis
}
\examples{

data(tpod)
test=gwas2(y=y,gen=gen[,1:240],fam=fam,chr=chr[1:12])
par(mfrow=c(2,1))

# Example Manhattan 1
SIGNIF = 1+(2*test$PolyTest$lrt>4.9)
plot(x=test,pch=SIGNIF+3,lwd=SIGNIF,main="Example 2")

# Example Manhattan 2
plot(x=test,main="Example 3",pch=20,lwd=2)
Kern = ksmooth(1:240,test$PolyTest$lrt,bandwidth=1)
lines(Kern,type="l",lwd=2)

}
\keyword{ NAM }
\keyword{ GWAS }
