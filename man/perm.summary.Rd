\name{perm.summary}
\alias{perm.summary}
\title{
This function computes various quantiles of the Null distributions of test statistics derived via permutation
}
\description{
Computes quantiles of null distributions. Includes the option of storing  values in tails which is useful
for combining results split over several processors.
}
\usage{
perm.summary(obj, sigma.inv, pval, nval = NULL, double.sided = TRUE, tails = FALSE)
}
\arguments{
  \item{obj}{
List - a perms object produced by the function perms.
}
  \item{sigma.inv}{
Sparse symmetricMatrix (dsCMatrix) - fitted sparse inverse covariance matrix.
}
  \item{pval}{
scalar - the pval to be used in the tests.
}
  \item{nval}{
scalar - the number of values in the tails of the null distribution to be stored. Should be >= trunc(nperm*pval/2+0.5)+2 .
}
  \item{double.sided}{
Logical - TRUE for double sided tests FALSE otherwise.
}
  \item{tails}{
logical - TRUE to store values in tails.
}
}
\details{
Null distribution summaries are produced for each contrast and the statistics for testing differential expression, differential connection and
the components of the overall test T. The actual nperm values derived by permutation are stored for the overall test T. 
}
\value{
A list with components
\item{deq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the (differential expression) contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{dcq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the neigbhour corrected (differential connection) contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{Tcq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the individula components of the overal test statistic T for the contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{Tq }{A list of length equal to the number of contrasts. Each component of this list is a vector of length equal to the total number of permutations and
consists of the realisations of the statistic T unde the Null hypothesis.}
\item{nperm }{ scalar - the total number of permutations.}
\item{seed }{ scalar - the seed used to create the permutation distributions.}
\item{tails}{tails - logical denoting whether obsevations in the tails are stored in the perm.summary object.}
\item{pval}{the pvalue requested.}
\item{ctb0}{ p by s matrix of observed contrast values, where s is the number of contrasts. }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Harri Kiiveri
}
\note{
If the total number of permutations have been done on one processor and you are not wanting to 
combine results from several processors use tails=FALSE.  For use when combining results from seversl processors, tails=TRUE and nval must be computed 
using nperm equal to the Total number of permutaions across all processors.  Note that nperm*pval/2 
should not be less than one, otherwise more permutations will be required to get sufficient observations in the tails.
}
\seealso{
\code{\link{combine.quantiles}}
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
# get summaries of permutation distributions
tmp<-perm.summary(res,sigma.inv=sinv,pval=.05,nval=NULL,double.sided=TRUE,tails=FALSE)
# quantiles of first contrast
tmp$deq[[1]]
# plot observed contrast 1 values against critical values
plot(tmp$ctb0[,1],ty="l",ylim=,xlab="variable",ylab="contrast value")
lines(tmp$deq[[1]][1,],col=2,lty=3)
lines(tmp$deq[[1]][5,],col=2,lty=3)
points(tmp$ctb0[,1])
}
\keyword{multivariate }
\keyword{ distribution }
