significance.tests <-
function(statistics.object,perm.summary.object,type=c("de","dc","cT","T"),top=20,gn)
  {
  if(perm.summary.object$tails)return("Error: perm.summary.object must have tails=FALSE","\n")
  pval<-perm.summary.object$pval
  results<-list()
  if(type[1]=="T"){
  a<-perm.summary.object$Tq
  s<-length(a)
  for(i in 1:s){
  cat("Contrast ",i,"\n")
  b<-quantile(a[[i]],probs=c(0.95,0.99))
  b<-c(statistics.object$T[i],b)
  names(b)<-c("overall test T",".05 level",".01 level")
  results[[i]]<-b
  print(b)
  cat("","\n")
  }
  return(results)
  }
  else{
  if(type[1]=="de"){s0<-statistics.object$G
  a<-perm.summary.object$deq}
  if(type[1]=="dc"){s0<-statistics.object$GN
  a<-perm.summary.object$dcq}
  if(type[1]=="cT"){s0<-statistics.object$cT
  a<-perm.summary.object$Tcq}
  s<-ncol(s0)
  #
  for(i in 1:s){
  cat("Contrast ",i,"\n")
  results[[i]]<-top.table(s0[,i],a[[i]],type=type[1],gn,top=top,pval)
  cat("","\n")
  }
  return(results)
  }
  }

