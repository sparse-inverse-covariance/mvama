qfunc <-
function(x,nval)
    {
    p<-length(x)
    i<-trunc(c(0.25,.5,.75)*p)
    x<-sort(x)
    res<-c(x[1:nval],x[i],x[(p-nval+1):p])
    res
    }

