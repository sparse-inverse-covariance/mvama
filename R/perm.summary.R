perm.summary <-
function(obj,sigma.inv,pval,nval=NULL,double.sided=TRUE,tails=FALSE)
{
# computes quantiles of null distributions and stores nval values in tails used
# for combining results split over several processors
if(class(obj)!="perms")return("Error: obj must be of class perms")
if(tails&is.null(nval))return("Error: nval must be specified")
s<-length(obj[[1]])
deq<-list()
dcq<-list()
Tcq<-list()
Tq<-list()
#
if(!tails){
for(i in 1:s){
P<-obj$CTB[[i]]
# quantiles for differential expression
deq[[i]]<-get.quantiles(P,pval,double.sided)
# get quantiles for differential connection
P<-sigma.inv%*%P
dcq[[i]]<-get.quantiles(P,pval,double.sided)
# get quantiles for components of T
P<-obj$CTB[[i]]*P
Tcq[[i]]<-get.quantiles(P,pval,double.sided)
# get all realisations for T
T<-colSums(P)
T<-matrix(T,nrow=1)
#Tq[[i]]<-get.quantiles(T,pval,double.sided=FALSE)
Tq[[i]]<-T
}
res<-list(deq=deq,dcq=dcq,Tcq=Tcq,Tq=Tq,nperm=obj$perm,seed=obj$seed,tails=tails,pval=pval,ctb0=obj$ctb0)
class(res)<-c("perm.summary")
return(res)
}
else{
for(i in 1:s){
P<-obj$CTB[[i]]
# tails for differential expression
deq[[i]]<-apply(P,1,qfunc,nval)
# get tails for differential connection
P<-sigma.inv%*%P
dcq[[i]]<-apply(P,1,qfunc,nval)
# get tails for components of T
P<-obj$CTB[[i]]*P
Tcq[[i]]<-apply(P,1,qfunc,nval)
# get tails for T
T<-colSums(P)
T<-matrix(T,nrow=1)
#Tq[[i]]<-apply(T,1,qfunc,nval)
Tq[[i]]<-T
}
}
res<-list(deq=deq,dcq=dcq,Tcq=Tcq,Tq=Tq,nperm=obj$nperm,seed=obj$seed,tails=tails,pval=pval,ctb0=obj$ctb0)
class(res)<-c("perm.summary")
return(res)
}

