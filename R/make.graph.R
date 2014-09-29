make.graph <-
function(sinv,gn=NULL,ide=NULL,idc=NULL,de="red",dc="green",both="yellow",neither="grey")
{
# can also use colour numbers - for a list of colors known to R type colors()
cols<-as.character(c(de,dc,both,neither))
adj<-sinv
if(!is.null(gn))colnames(adj)<-gn
adj@x<-rep(1,length(adj@x))
if(!is.null(gn))colnames(adj)<-gn
if(is.null(ide))ide<-rep(0,ncol(sinv))
if(is.null(idc))idc<-rep(0,ncol(sinv))
colr<-ide+2*idc   # nothing - white, differential expression= black, differential connection=red, both=green
# change colours
colr[colr==1]<-cols[1]   # differential expression red
colr[colr==2]<-cols[2]   # differential connection green
colr[colr==3]<-cols[3]   # both de and dc - yellow
colr[colr==0]<-cols[4]   # neither - grey
g<-graph.adjacency(adj,diag=FALSE,mode="undirected",add.colnames="label")
V(g)$color<-colr
g
}

