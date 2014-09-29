get.paths <-
function(g,from,to)
{
# wrapper to deal with character labels in from and to and indexing
if(is.null(V(g)$label))return("Error: g must be labelled")
gn<-V(g)$label
if(is.character(from)){
i<-match(from,gn)
if(is.na(i))return("Error: from label invalid")
}
else {i<-from }
if(is.character(to)){
j<-match(to,gn)
if(any(is.na(j)))return("Error: to label invalid")
}
else{j<-to }
a<-get.all.shortest.paths(g,from=i[1]-1,to=j-1)
if(length(a)==0){
a<-list()
return(a)}
for(i in 1:length(a)){
a[[i]]<-a[[i]]+1
a[[i]]<-V(g)$label[a[[i]]]
}
a
}

