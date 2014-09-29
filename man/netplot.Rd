\name{netplot}
\alias{netplot}
\title{
Plot an undirected graph represented as an igraph object
}
\description{
Plot an undirected graph represented as an igraph object
}
\usage{
netplot(g, layout = layout.kamada.kawai, seed = 123, ...)
}
\arguments{
  \item{g}{
an igraph object
}
  \item{layout}{ see the igraph documention for layout options

}
  \item{seed}{
random seed to intialise layout algorithm.
}
  \item{\dots}{
other options for plotting an igraph object. See the example.
}
}
\details{
A convenience function for plotting the igraph objects used in this package.
A wrapper for plot(..)
}
\value{
Silently plots to the chosen device.
}
\references{
See the igraph reference manual.
}
\author{
Harri Kiiveri
}
\note{
Some useful arguments are vertex.label.color and edge.color
}
\seealso{
\code{\link{make.graph}}
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
# other options
netplot(g,edge.color="red")
netplot(g,edge.color="red",vertex.label.color="black")
}
}
\keyword{ graphs }
\keyword{hplot}
