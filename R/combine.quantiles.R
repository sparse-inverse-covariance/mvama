combine.quantiles <-
function(obj.list)
  {
  # computes p values and quantiles by combining from a list of perm.summary objects
  n<-length(obj.list)
  w<-rep(0,n)
  test<-0
  seed<-w
  for(i in 1:n){
  w[i]<-obj.list[[i]]$nperm
  seed[i]<-obj.list[[i]]$seed
  test<-test+as.numeric((class(obj.list[[i]])[1]=="perm.summary")&(obj.list[[i]]$tails))
  }
  if(test<n)return("Error:objects must be class perm.summary with tails=TRUE")
  # compute total number of permutations
  nsim<-sum(w)
  w<-w/nsim
  # set up needed parameters
  d<-dim(obj.list[[1]][[1]][[1]])
  s<-length(obj.list[[1]][[1]])    # no of contrasts
  nr<-d[1]
  nc<-d[2]
  nval<-(nr-3)/2
  pval<-obj.list[[1]]$pval/2      #  two sided tests
  #
  il<-1:nval
  im<-(nval+1):(nval+3)
  iu<-(nval+4):nr
  # loop over statistics and obj.list components
  stor<-list()
  alist<-list()
  for(i in 1:3){
  #if(i==4){nc<-1
  # pval<-pval*2}  #  two sided tests except for T
  for(k in 1:s){
  # cat(i,k,"\n")
  M<-matrix(0,nrow=3,ncol=nc)
  L<-NULL
  U<-L
  for(j in 1:n){
  Q<-obj.list[[j]][[i]][[k]]
  M<-M+w[j]*Q[im,,drop=FALSE]
  L<-rbind(L,Q[il,,drop=FALSE])
  U<-rbind(U,Q[iu,,drop=FALSE])
  }
  LE<-apply(L,2,sort)
  UE<-apply(U,2,sort)
  nu<-nrow(UE)
  LE<-LE[il,,drop=FALSE]
  UE<-UE[(nu-nval+1):nu,,drop=FALSE]
  # get pvals
  resl<-getjg(nsim,pval)
  resu<-getjg(nsim,1-pval)
  resu$j<-resu$j-(nsim-nval)
  # lower quantiles
  lq<-(1-resl$g)*LE[resl$j,]+resl$g*LE[resl$j+1,]
  # upper quantiles
  uq<-(1-resu$g)*UE[resu$j,]+resu$g*UE[resu$j+1,]
  #
  tmp<-rbind(lq,M,uq)
  rownames(tmp)<-c(pval,.25,.5,.75,1-pval)
  alist[[k]]<-tmp
  }
  stor[[i]]<-alist
  }
  # do T statistic separately - save all realisations
  for (k in 1:s){
  a<-NULL
  for(j in 1:n){
  a<-cbind(a,obj.list[[j]][[4]][[k]])
  }
  alist[[k]]<-a
  }
  res<-list(deq=stor[[1]],dcq=stor[[2]],Tcq=stor[[3]],Tq=alist,nperm=nsim,seed=seed,tails=FALSE,pval=2*pval,ctb0=obj.list[[1]]$ctb0)
  class(res)<-c("perm.summary", "combined")   # same as perm.summary object except seed is a vector of seeds
  return(res)
  }

