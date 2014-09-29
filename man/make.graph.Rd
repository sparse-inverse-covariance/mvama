\name{make.graph}
\alias{make.graph}
\title{
create an igraph object
}
\description{
This is a convenience function to facilitate the creation of graphical objects 
useful in displaying results in the context of a gene network
}
\usage{
make.graph(sinv, gn = NULL, ide = NULL, idc = NULL, de = "red", dc = "green", both = "yellow", neither = "grey")
}
\arguments{
  \item{sinv}{
sparse symmetric matrix. Can be the estimated inverse covariance matrix or a symmetric adjacency matrix.
}
  \item{gn}{
real vector - gene/variable names to place on vertices.
}
  \item{ide}{
real vector - (0,1) indicator vector identifying differentially expressed genes.
}
  \item{idc}{
real vector - (0,1) indicator vector identifying differentially connected genes.
}
  \item{de}{
character or integer - vertex colour for differentially expressed genes.
}
  \item{dc}{
character or integer - vertex colour for differentially connected genes.
}
  \item{both}{
character or integer - vertex colour for differentially connected genes.
}
  \item{neither}{
character or integer - vertex colour for genes which are neither differentially expressed or differentially connected.
}
}
\details{
Generates an igraph object corresponding to an undirected graph with coloured and labelled  vertices.
}
\value{
Returns a labelled and coloured igraph object
}
\references{
See the igraph package reference manual.
}
\author{
Harri Kiiveri
}
\note{
If gn is null vertex labels are taken to be the column names of sinv. Colours mey be defined via text
(for a list of colors known to R type colors() in R) or can also be defined by using integers. 
}
\seealso{
\code{\link{netplot}}
}
\examples{
\dontrun{
library(Matrix)
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# make up stuff about differential expression and differential connection
de<-c(0,1,0,1,1)
dc<-c(0,0,0,0,1)
# create graph
g<-make.graph(sinv,gn=c("a","b","c","d","e"),ide=de,idc=dc)
# plot
netplot(g,vertex.label.color="black")
}
}
\keyword{ graphs }
