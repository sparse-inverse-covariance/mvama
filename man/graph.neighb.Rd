\name{graph.neighb}
\alias{graph.neighb}
\title{
Produces neighbourhoods of vertices of a graph 
}
\description{
For each vertex in a specified set produces a subgraph coresponding to the 
vertices within a specified distance.
}
\usage{
graph.neighb(g, order = 3, nodes)
}
\arguments{
  \item{g}{
an igraph object produced by make.graph.
}
  \item{order}{
scalar - integer size of the neighbourhood.
}
  \item{nodes}{
vector - integers or characters defining the vertices for which to calculate neighbourhoods.
}
}
\details{
A wrapper for the igraph function graph.neighborhood.
}
\value{
A list of graphs.
}
\references{
See igraph reference manual.
}
\author{
Harri Kiiveri
}
\note{
Labels and colours are inherited from g.
}
\seealso{
\code{\link{get.paths}}
}
\examples{
\dontrun{
library(Matrix)
#generate sigma.inv
sinv<-diag(5)/2
sinv[1,2]=sinv[2,3]=sinv[3,4]=sinv[4,5]=sinv[3,5]=0.1
sinv<-sinv+t(sinv)
sinv<-as(sinv,"dsCMatrix")
# make up stuff about differential expression and differential connection
de<-c(0,1,0,1,1)
dc<-c(0,0,0,0,1)
# create graph
g<-make.graph(sinv,gn=c("a","b","c","d","e"),ide=de,idc=dc)
netplot(g)
res<-graph.neighb(g,order=1,"c")
# plot neighbourhood
netplot(res[[1]])
}
}
\keyword{ graphs }

