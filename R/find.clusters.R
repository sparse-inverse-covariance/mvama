find.clusters <-
    function(g,nodes )
{
    if(is.character(nodes)){
        gn<-V(g)$label
        nodes<-match(nodes,gn)
    }
    if(any(is.na(nodes)))return("Error: invalid nodes")
    # find clusters
    f<-subgraph(g,nodes)  # vertex numbering now starts from 1 in igraph
    res<-clusters(f)
    #create subgraphs of each cluster
    clust<-list()
    for(i in 1:res$no){
        clust[[i]]<-subgraph(f,(which(res[[1]]==i)))
    }
    list(gs=f,clust.list=clust,csize=res$csize,cno=res$no)
}
