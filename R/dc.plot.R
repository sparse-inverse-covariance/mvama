dc.plot <-
function(dc.gene.id,pdata,B,ctb0,...)
{
# produces a neighbour plot for differentially connected genes
# dc.gene.id is id of a differentially connected gene
# pdata is a p by nsim matrix of contrast realisations from the null distribution
# B is the matrix of regression coefficients (diagonal=variances) - !!with column names
# ctb0 is the vector of observed contrasts
i<-which(B[dc.gene.id,]!=0)
i0<-which(i==dc.gene.id)
ii<-c(i[i0],i[-i0])
# point for predicted contrast value based on neighbours
b<-B[dc.gene.id,]
b[dc.gene.id]<-0
ghat<-sum(b*ctb0)
#
ctb0<-ctb0[ii]
tmp<-t(pdata[ii,])
colnames(tmp)<-colnames(B)[ii]
# put response gene first
yhigh<-apply(tmp,2,max)
ylow<-apply(tmp,2,min)
yhigh<-max(yhigh,ctb0)
ylow<-min(ylow,ctb0)
del<-(yhigh-ylow)/10
boxplot(tmp,ylim=c(ylow-del,yhigh+del),ylab="contrast value",...)
lines(c(1.5,1.5),c(ylow-del,yhigh+del),lty=3,col=1)
lines(ctb0,col=2)
points(ctb0,col=2,pch=19)
nc<-ncol(tmp)
x<-1:nc
y1<-rep(yhigh+del/2,nc)
y2<-rep(ylow-del/2,nc)
# plot regression coeffs and sd of dc.gene.id
text(x,y1,colnames(tmp))
# plot gene names
text(x,y2,round(B[dc.gene.id,ii],2))    # changed from i to ii  ------
# plot point for predicted contrast value based on neighbours
points(1,ghat,pch=19,col=3)
}

