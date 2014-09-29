\name{sigma.inv.approx}
\alias{sigma.inv.approx}
\title{
Fast approximation to sigma inverse
}
\description{
Uses regression to compute a fast approximation to the inverse covariance matrix
}
\usage{
sigma.inv.approx(R, a, index = 1:nrow(a), symmetric = FALSE)
}
\arguments{
  \item{R}{
n by p matrix - mean corrected residual data matrix
}
  \item{a}{
sparse symmeric matrix - either a sparse symmetric adjacency matrix defining the pattern of 
non zeroes in sigma inverse or a fitted inverse covariance matrix.
}
  \item{index}{
vector - specifying which rows of the approximation to compute.
}
  \item{symmetric}{
logical - if TRUE return a symmetric matrix.
}
}
\details{
For each row k in index, a linear regression of variable k on its neighbours 
defined by the adjacency matrix is computed. The regression coefficients and 
residual variance are then transformed into the non zero entries of the kth row
of the approximate sigma inverse.
}
\value{
A sparse matrix of class dgCMatrix (or dsCMatrix if symmetric=TRUE). 
}
\references{

}
\author{
Harri Kiiveri
}
\note{
If symmetric=TRUE the resulting esitmate is (a+t(a))/2 where a is the estimate for symmetric=FALSE.
Note that the estimate may not be positive definite.
}

\seealso{

}
\examples{
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
colnames(sinv)<-c("a","b","c","d","e")
# generate data
res<-Simulate(sinv,nsamp=10000)
X<-res$data[[1]]
# mean correct
R<-scale(X,center=TRUE,sc=FALSE)
tmp<-sigma.inv.approx(R,sinv,symmetric=TRUE)
# check results
sinv # true values
tmp  # approximation
}
\keyword{ multivariate }
\keyword{ models }
