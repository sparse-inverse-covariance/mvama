get.statistics <-
function(obj,sigma.inv)
            {
            # compute statistics from a perms or perm.summary  object
            G<-obj$ctb0
            GN<-sigma.inv%*%G
            cT<-G*GN
            T<-colSums(cT)
            list(G=G,GN=GN,cT=cT,T=T)
            }

