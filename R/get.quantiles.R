get.quantiles <-
function(P,pval,double.sided)
{
# copmute quantiles of each row of p by nperm matrix P
if( double.sided)pval<-pval/2
a<-apply(P,1,quantile,probs=c(pval,.25,0.5,.75,1-pval))
a
}

