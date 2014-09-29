\name{get.statistics}
\alias{get.statistics}
\title{
Computes observed statistics for significance testing
}
\description{
Computes the observed quantities required to test for differential expression and 
differential connection both overall and componentwise.
}
\usage{
get.statistics(obj, sigma.inv)
}
\arguments{
  \item{obj}{
a "perms" object or a "perm.summary" object.
}
  \item{sigma.inv}{
sparse symmetric matrix - the fitted inverse covariance matrix.
}
}
\details{
If g denotes a p by 1 vector of observed contrasts then the neighbour corrected contrasts are defined by
gn=sigma.inv\%*\%g. The components of T are the components g[i]*gn[i]. The statistic T is t(g)\%*\%sigma.inv\%*\%g.
}
\value{
A list with components
\item{G }{a p by s matrix of observed contrast values, where s is the number of contrasts.}
\item{GN }{a p by s matrix of neighbour corrected contrast values.}
\item{cT }{a p by s matrix of the observed components of the T statistic for each contrast.}
\item{T }{a vector of length s giving the T statistic for each contrast.}
}
\references{
Kiiveri(2011) Multivariate analysis of microarray data: differential expression and differential connection. BMC Bioinformatics.
}
\author{
Harri Kiiveri
}
\note{
If experimenting with different estimates of sigma inverse, note that this function will accept sigma.inv which is 
not positive definite or symmetric.
}
\seealso{
\code{\link{significance.tests}}
}
\examples{
library("Matrix")
# small 3 group problem with 5 variables and 150 samples
# generate mean matrix
a<-as.numeric(gl(3,50,150))
b<-c(0,1,0,1,1)
M<-outer(a,b,"*") # variables 2,4 and 5 have nonzero and unequal group means
# design matrix
D<-matrix(0,nrow=150,ncol=3)
D[1:50,1]<-1
D[51:100,2]<-1
D[101:150,3]<-1
# generate data
X<-matrix(rnorm(5*150),ncol=5)+M
# contrast matrix for testing equality of group2 and 1 and group 3 and group 1 means
C<-matrix(c(1,-1,0,1,0,-1)/2^0.5,ncol=2)
# generate permutation distributions
res<-perms(C,X,D,nperm=100,seed=123)
#------------------------------------------------------------------------------
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
tmp<-get.statistics(res,sinv)
tmp
}
\keyword{ htest }
