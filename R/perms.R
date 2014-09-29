perms <-
function(C,X,D,nperm=500,seed=123)
   {
   # computes permutation test for the contrast C'B - C is r by s B is r by p
   # X is n by P
   # D is n by r
   # reparameterise
   # compute observed contrasts 
   p<-ncol(X)
   res<-lm.fit(D,X)$coefficients
   C.sav<-C
   if(!is.matrix(C))C<-matrix(C,ncol=1)
   r<-nrow(C)
   s<-ncol(C)
   CB<-t(C)%*%res
   set.seed(seed)
   n<-nrow(X)
   p<-ncol(X)
   r<-ncol(D)
   s<-ncol(C)
   #C<-diag(r)-C%*%t(C)        
   C<-diag(r)- C%*%solve(crossprod(C))%*%t(C)
   Q<-svd(C)$u
   # avoid sign changes in contrasts
   Q[,(r-s+1):r]<-C.sav
   # compute new design matrix
   # D<-D%*%Q                  
   D<-D%*%t(solve(Q))
   # fitted values under null C'B=0
   res<-lm(X~D[,1:(r-s)]-1)
   FV<-res$fitted.values
   R<-res$residuals
   qr<-qr(D)
   # do permutations
   CTB<-list()
   rms<-r-s
   for(k in 1:s){
   CTB[[k]]<-matrix(0,nrow=p,ncol=nperm)
   }
   for(i in 1:nperm){
   if(i%%500==0)cat("completed ",i,"\n")
   ii<-sample(1:n)
   XP<-FV+R[ii,]
   for(k in 1:s){
   CTB[[k]][,i]<-qr.coef(qr,XP)[rms+k,]
   }      
   }
   # return list of contrast matrix values under null and new design matrix
   # required to compute tests (using last r-s values)
   # DQ is transformed design matrix, ctb0 is p by s observed values of contrast(s) for the data set
   res<-list(CTB=CTB,DQ=D,ctb0=t(CB),nperm=nperm,seed=seed)
   class(res)<-c("perms")
   res
   }

