\name{prediction.plot}
\alias{prediction.plot}
\title{
Plots the contrast values of a variable and its neighbours
}
\description{
Plots the contrast values of a variable and its neighbours in the context of the null distributions of each of these contrasts.
}
\usage{
prediction.plot(gene.id, sinv, obj, contrast = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.id}{
character or scalar - identifier of the variable of interest.
}
  \item{sinv}{
sparse symmetric matrix - the fitted inverse covariance matrix.
}
  \item{obj}{
a "perms" object or a "perm.summary" object produced by the function perms or perm.summary.
}
  \item{contrast}{
scalar - identifying which contrast to use in plot.
}
}
\details{
The observed contrast values of the variable and its neighbours are plotted on 
top of boxplots of their null distributions. The observed values are plotted in red.
The predicted value based on the weighted linear combination of neighbour values is 
ploted in green. Regression coefficients are also included at the bottom of the plot.
Predictor variables are to the right of the dotted vertical line. The number on the 
bottom left is the residual variance in the fitted regression model.
}
\value{
Produces a plot.

}
\references{

}
\author{
Harrri Kiiveri
}
\note{
Variable labels are taken from the column names of sinv and are also plotted.
}
\seealso{
\code{\link{graph.neighb}}
}
\examples{
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
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
colnames(sinv)<-c("a","b","c","d","e")
# generate permutation distributions
res<-perms(C,X,D,nperm=100,seed=123)
prediction.plot(gene.id="c",sinv,obj=res,contrast=1)
}
\keyword{ hplot }
