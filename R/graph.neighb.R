graph.neighb <-
function(g,order=3,nodes)
{
# for each element of nodes  produces subgraph defined by all connected vertices
# within a distance order of each vertex in nodes
if(is.character(nodes)){
gn<-V(g)$label
nodes<-match(nodes,gn)
}
if(any(is.na(nodes)))return("Error: invalid nodes")
nodeids<-nodes
tmp<-graph.neighborhood(g,order=order,nodes=nodeids-1) # NOTE usual vertex label-1
tmp
}

