\name{combine.quantiles}
\alias{combine.quantiles}
\title{
Combines quantiles obtained from several independent sets of permutations
}
\description{
Gets  p values and quantiles for the combined set of permutations by combining from a list of 
perm.summary objects for the component pieces. This function is useful for combining results
obtained by splitting the permuations over different processors or runs to save memory and time.
}
\usage{
combine.quantiles(obj.list)
}
\arguments{
  \item{obj.list}{ -  a list of perm.summary objects produce by the function perm.summary 
with tails=TRUE.
}
}
\details{
Estimates the 25\% 50\% and 75\% quantiles for the combined set of permutations
by averaging the values for each independent component set. The 0.5*pval*100\% and
(1-0.5pval)*100\% quantiles are calculated exactly for the combined set given the
 observations in the tails of each component data set.
}
\value{
A list with components
\item{deq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the (differential expression)  contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{dcq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the neigbhour corrected (differential connection) contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{Tcq }{A list of length equal to the number of contrasts. Each component of this list is a 5 by p
matrix  of critical values for testing the significance of the individula components of the overal test statistic T for the contrast. The rows correspond to p=c(pval/2,0.25,0.5,0.75,1-pval/2) and the columns to the p
variables or genes.}
\item{Tq }{A list of length equal to the number of contrasts. Each component of this list is a vector of length equal to the total number of permutations and
consists of the realisations of the statistic T unde the Null hypothesis.}
\item{nperm }{ scalar - the total number of permutations}
\item{seed }{ scalar or vector - the seed(s) used to create the permutation distributions}
\item{tails}{tails - logical denoting whether obsevations in the tails are stored in the perm.summary objects}
\item{pval}{the pvalue requested}
\item{ctb0}{ p by s matrix of observed contrast values, where s is the number of contrasts. }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Harri Kiiveri
}
\note{
The perm.summary objects in obj.list must be created with common values of sigma.inv,pval,nval,double.sided.
A perm.summary object is produced with a vector of values for seed.
}
\seealso{
\code{\link{perm.summary}}
}
\examples{
library(Matrix)
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
# make up sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# generate data
res<-Simulate(sigma.inv=sinv,M=M,nsamp=150,nsim=1,seed=123)
X<-res$data[[1]]
# contrast matrix for testing equality of group2 and 1 and group 3 and group 1 means
C<-matrix(c(1,-1,0,1,0,-1)/2^0.5,ncol=2)
#
#-------------------------------------------------------------------------------
# generate permutation distributions - get quantiles for 200 permutations by combining two sets of
# 100 permutations each
# on processor 1 say
nval<-trunc(200*.05/2+0.5)+2
res1<-perms(C,X,D,nperm=100,seed=123)
a1<-perm.summary(res1,sigma.inv=sinv,pval=.05,nval=nval,double.sided=TRUE,tails=TRUE)
# on processor 2
nval<-trunc(200*.05/2+0.5)+2
res2<-perms(C,X,D,nperm=100,seed=456)
a2<-perm.summary(res2,sigma.inv=sinv,pval=.05,nval=nval,double.sided=TRUE,tails=TRUE)
# make list
a<-list()
a[[1]]<-a1
a[[2]]<-a2
# combine quantiles
tmp<-combine.quantiles(a)
tmp$deq[[1]]
tmp$nperm

#------------------------------------------------------------------------------
}
\keyword{ multivariate }
\keyword{ distribution}
