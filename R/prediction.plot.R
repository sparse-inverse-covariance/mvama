prediction.plot <-
function(gene.id,sinv,obj,contrast=1)
{
s<-contrast
# compute matrix of regression coefficients
if((class(obj)[1]!="perms")&(class(obj)[1]!="perm.summary"))return("Error: obj must be class perms or perm.sumary" ,"\n")
gn<-colnames(sinv)
d<-diag(sinv)
d<-Diagonal(n=ncol(sinv),1/d)
B<-d%*%sinv
B<--B
diag(B)<-diag(d)
# gene.id can be name or integer - check
if(is.character(gene.id))rr<-which(gn==gene.id)
else rr<-gene.id
# extract required components
if(class(obj)[1]=="perms"){
pdata<-obj$CTB[[s]]
ctb0<-obj$ctb0[,s]
dc.plot(rr,pdata,B,ctb0)
}else{
# surprise - it works here as well - although no outliers available to be plotted
a<-obj$deq[[s]]
ctb0<-obj$ctb0[,s]
aa<-a
del<-aa[4,]-aa[2,]
aa[1,]<-aa[2,]-1.5*del  # make whiskers extend to 1,5 times interquartile range beyond box
aa[5,]<-aa[4,]+1.5*del
dc.plot(rr,t(aa),B,ctb0,range=0)
}
}

