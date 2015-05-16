sigma.inv.approx <-function(R,a,index=1:nrow(a),symmetric=FALSE){
    n<-nrow(R)
    sigma.inv<-a
    sigma.inv<-as(sigma.inv,"dgCMatrix")
    for(i in index){
        k<-which(a[i,]!=0)
        if(length(k)>1){
            r<-length(k)
            v<-rep(0,r)
            i0<-which(k==i)
            tmp<-R[,k]
            res<-lm(tmp[,i0]~tmp[,-i0]-1)
            b<-res$coefficients
            s<-sum(res$residuals^2)/n
            v[i0]<-1/s
            v[-i0]<- -b/s
            sigma.inv[i,k]<-v
        }
        else{
            sigma.inv[i,i]<-1/(sum(R[,i]^2)/n)
        }
    }
    if(symmetric){sigma.inv<-(sigma.inv+t(sigma.inv))/2   # change class to symmetric 
                  sigma.inv <-forceSymmetric(sigma.inv) 
              }
    sigma.inv
}
