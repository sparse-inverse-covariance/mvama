top.table <-
function(stats,Q,type,gn,top=20,pval)
  {
  # get median and robust standard deviation
  m<-Q[3,]
  sd<-(Q[4,]-Q[2,])/1.3489
  # scale observed test statistics by median and interquartile range
  sctb0<-(stats-m)/sd
  # indicator for differential expression
  ide<-(stats<Q[1,])| (stats>Q[5,])
  nde<-sum(ide)
  #  print out table of differentially expressed genes
  k<-which(ide)
  datde<-data.frame(gn[k],stats[k],Q[1,k],Q[5,k],sctb0[k],sd[k])
  colnames(datde)<-c("gene","statistic","lower quantile","upper quantile", "scaled statistic","robust sd")
  rownames(datde)<-NULL
  kk<-order(abs(datde[,5]),decreasing=TRUE) # order table by absolute scaled statistic
  #print top ones
  if(type=="de")cat("Top",top," differentially expressed genes - p =",pval,"\n")
  if(type=="dc")cat("Top",top," differentially connected genes - p =",pval,"\n")
  if(type=="cT")cat("Top",top," genes with large T statistic components - p =",pval,"\n")
  print(datde[kk[1:top],])
  cat("Total number significant",nde,"\n")
  cat("","\n")
  #ids of top genes
  deids<-k[kk[1:top]]
  list(topids=deids,allids=k,top.table=datde[kk[1:top],],nsig=nde)
  }

