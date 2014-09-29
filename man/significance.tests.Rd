\name{significance.tests}
\alias{significance.tests}
\title{
This function calculates significance tests for contrasts
}
\description{
Calculates significance tests for contrasts, neighbour corrected contrasts, components of the statistic T and for T.
}
\usage{
significance.tests(statistics.object, perm.summary.object, type = c("de", "dc", "cT", "T"), top = 20, gn)
}
\arguments{
  \item{statistics.object}{
list object - produced by get.statistics
}
  \item{perm.summary.object}{
list object - produced by perm.summary
}
  \item{type}{
character - the type of test to compute
}
  \item{top}{
scalar - output the top most significant variables/genes
}
  \item{gn}{
character vector - names of the variables/genes (length(p))
}
}
\details{
For type = "de","dc" or "cT" writes a table to the screen for each contrast, see below. For type="T", for each contrast writes out the observed T value
and the .05 and .01 levels of the upper tail of the  Null distribution.

}
\value{
A list with s components where s is the number of contrasts . For type= "de","dc" or "cT",each of these components is itself a list with components
\item{ topids}{id numbers of the top significant genes}
\item{allids }{id numbers of all significant genes}
\item{top.table }{top by 6 table with columns gene name, statistic, lower quantile, upper quantile, scaled statistic, robust standard deviation}
\item{nsig }{total number of significant genes}
For type="T" each list component is a vector of length 3 with elements, observed T test value,
 .05 level of Null distribution and .01 level of Null distribution.
}
\references{
Kiiveri, H. T.(2010) Multivariate analysis of microarray data: differential expression and differential connection.
}
\author{
Harri Kiiveri
}
\note{
The lower and upper quantiles correspond to  pval/2 and (1-pval/2) where pval is obtained from the perm.summary object. The robust standard
deviation is the interquartile range divided by 1.3489. The statistic is median corrected and scaled by the robust standard deviation.
Ranking is on the basis of the scaled statistic. The should be OK for type="de", but is not necessarily ideal for type="dc" or "cT". Future versions
may estimate the tail behaviour of the null distribution and use the tail pvalues to rank the statistics. The test for T is a one tailed test.
Other p values for this test can be obtained by computing the quantiles T from the realisations stored in a "perms" or a "perm.summary" object.
}
\seealso{
\code{\link{perms}},\code{\link{perm.summary}},\code{\link{get.statistics}}
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
#generate sigma.inv - illustration only - not consistent with data
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# get summaries of permutation distributions
tmp<-perm.summary(res,sigma.inv=sinv,pval=.05,nval=NULL,double.sided=TRUE,tails=FALSE)
# compute statistics
ress<-get.statistics(tmp,sinv)
# do significance tests
# overall test T
rT<-significance.tests(ress,tmp,type=c("T"),top=3,gn=c("a","b","c","d","e"))
# components of T
rcT<-significance.tests(ress,tmp,type=c("cT"),top=3,gn=c("a","b","c","d","e"))
# plot components of T to identify possible interesting connections for contrast 1
q<-tmp$Tcq[[1]]  # quantiles of NULL distribution of components of T
v<-ress$cT[,1]# observed components of T
# plot components of T and significance levels
plot(q[1,],ty="l",lty=3,col=2,ylim=c(-.1,1.1),ylab="cT value",xlab="variable")
lines(q[5,],col=2,lty=3)
lines(v)
points(v)
# differential expression
rde<-significance.tests(ress,tmp,type=c("de"),top=3,gn=c("a","b","c","d","e"))
# differential connection
rdc<-significance.tests(ress,tmp,type=c("dc"),top=3,gn=c("a","b","c","d","e"))
}
\keyword{ htest }
