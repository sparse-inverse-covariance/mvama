\name{Simulate}
\alias{Simulate}
\title{
Simulate matrix realisations from a multivariate normal distribution
}
\description{
Simulates realisations from a multivariate normal distribution with a covariance matrix with a specified (sparse) inverse.
}
\usage{
Simulate(sigma.inv, M = NULL, nsamp = 100, nsim = 1, seed = 123)
}
\arguments{
  \item{sigma.inv}{
sparse symmetricMatrix - inverse covariance matrix.
}
  \item{M}{
Matrix - mean value matrix (nsamp by ncol(sigma.inv) Defaults to 0.
}
  \item{nsamp}{
scalar - number of rows of each data matrix (stored as nsamp by p matrix).
}
  \item{nsim}{
scalar - number of nsamp by p matrices to generate.
}
  \item{seed}{
scalar - seed for generating variates.
}
}
\details{
Row i of each generated data matrix is N(M[i,],sinv^-1). Uses sparse matrix calculations 
suitable for very high dimension p to do simulations.
}
\value{
A list with components
\item{data }{a list of length nsim, each component of which is a nsamp by p matrix.}
\item{n }{scalar - the number of  nonzero entries in the upper or lower triangle of sigma.inv.}
\item{naf }{scalar - the number of nonzeroes in the sparse cholesky factorisation of sigma.inv.}
\item{L }{ sparse matrix - the sparse cholesky facxtor of sigma.inv ordered to minimise fillin.}
\item{perm }{ vector - the permutation required to minimise fillin.}

}
\references{
}
\author{
Harri Kiiveri
}
\note{
Positive definiteness of sigma.inv is not checked.
}
\seealso{
}
\examples{
library("Matrix")
#generate sigma.inv 
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# generate data
res<-Simulate(sigma.inv=sinv,M=NULL,nsamp=10000,nsim=1,seed=123)
X<-res$data[[1]]
X<-scale(X,center=TRUE,sc=FALSE)
# check means
colMeans(X)
# check structure in inverse
solve(crossprod(X)/10000)
# true inverse
sinv
}
\keyword{ multivariate }
