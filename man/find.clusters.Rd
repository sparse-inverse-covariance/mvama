\name{find.clusters}
\alias{find.clusters}
\title{
Find clusters of connected vertices in a subgraph
}
\description{
This is a convenience function providing an interface to the igraph functions subgraph and clusters.
}
\usage{
find.clusters(g, nodes)
}
\arguments{
  \item{g}{
an igraph object - produced by make.graph
}
  \item{nodes}{
character or vector - a sequence of vertex numbers or vertex names to define the initial subgraph.
}
}
\details{
Finds connected components (clusters) of the subgraph of g defined by nodes and stores these components as igraph objects.
}
\value{
A list with components
\item{gs }{an igraph object - the subgraph of g defined by nodes.}
\item{clust.list }{a list with each component being an igraph object representing the graph of a cluster}
\item{csize }{ the vector of cluster sizes}
\item{cno }{the number of clusters}
}
\references{
See the igraph reference manual.
}
\author{
Harri Kiiveri
}
\note{
The clusters inherit the labels and colours of g.
}
\seealso{
\code{\link{make.graph}}
}
\examples{
\dontrun{
library(Matrix)
library(igraph)
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[4,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# make up stuff about differential expression and differential connection
de<-c(0,1,0,1,1)
dc<-c(0,0,0,0,1)
# create graph
g<-make.graph(sinv,gn=c("a","b","c","d","e"),ide=de,idc=dc)
res<-find.clusters(g,1:5)
#look at results
netplot(g)
netplot(res$gs)
netplot(res$clust.list[[1]])
netplot(res$clust.list[[2]])
res$csize
res$cno
}
}
\keyword{graphs}
