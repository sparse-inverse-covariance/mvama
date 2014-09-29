Simulate <-
function(sigma.inv,M=NULL,nsamp=100,nsim=1,seed=123)
{
p<-nrow(sigma.inv)
data<-list()
res<-Cholesky(sigma.inv,LDL=FALSE,perm=TRUE)
L<-expand(res)$L
P<-expand(res)$P
perm<-P@perm
iperm<-t(P)@perm
# permute a for minimum fillin
# a<-a[perm,perm]
L<-as(t(L),"dgCMatrix") # row rep = col rep of transpose
# simulate data with this covariance matrix in this order
set.seed(seed)
for(i in 1:nsim){
Y<-matrix(rnorm(p*nsamp),nrow=p)
X<-as.matrix(solve(L,Y))
# undo order
data[[i]]<-t(X[iperm,])
if(!is.null(M))data[[i]]<-data[[i]]+M
}
n<-length(sigma.inv@i);naf<-length(L@i)
list(data=data,n=n,naf=naf,L=L,perm=perm)
}

