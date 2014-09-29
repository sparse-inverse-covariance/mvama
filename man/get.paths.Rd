\name{get.paths}
\alias{get.paths}
\title{
Finds all shortest paths between a vertex and a set of vertices 
}
\description{
This is a convenience fumction which provides an interface to the igraph function get.all.shortest.paths.
}
\usage{
get.paths(g, from, to)
}
\arguments{
  \item{g}{
an igraph object
}
  \item{from}{
integer or character variable.
}
  \item{to}{
a vector of vertex numbers or labels.
}
}
\details{
See the igraph function get.all.shortest paths.
}
\value{
A list with each component being a sequence of vertices defining a path from the "from" vertex to a "to" vertex.
}
\references{
See the igraph reference manual.
}
\author{
Harri Kiiveri
}
\note{
}
\seealso{
\code{\link{make.graph}}, \code{\link{find.clusters}}
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
res<-get.paths(g,from="a",to=c("d","e"))
# print out paths
res
}
}
\keyword{ graphs }

