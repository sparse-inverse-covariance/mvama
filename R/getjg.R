getjg <-
function(n,pval)
    {
    m<-1-pval
    npm<-n*pval+m
    j<-floor(npm)
    g<-npm-j
    list(j=j,g=g)
    }

