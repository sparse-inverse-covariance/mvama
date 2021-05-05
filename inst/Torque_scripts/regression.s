# clean up Tuesday, March 29, 2011

if (!is.loaded("mpi_initialize")) { library("Rmpi") }

foldNumber<-NA
mpi.bcast.cmd(foldNumber <- mpi.comm.rank())

# Get command line input 
filename<-""
outdir<-""
restartdir<-""
restart<-FALSE
checkpointcycles<-100
exitoncheckpoint<-FALSE
kmax<-2
# get only input for this script
args=(commandArgs(TRUE))
# for each arguement evaluate the statement
# e.g. "filename=/"some_file.Rdata/""
for ( i in args ) {
  eval(parse(text=i))
}
# Default behaviour if certain arguements aren't specified
if( filename=="" ){
  mpi.close.Rslaves()
  mpi.quit()
}
if( outdir=="" ){
  outdir<-"./"
}
if( restart==TRUE && restartdir=="" ) {
  mpi.close.Rslaves()
  mpi.quit()
}

if( restart==TRUE ){
  print(paste("Restarting processing from checkpoint files in",restartdir,sep=" "))
}

myfunc<-function(){
  
  library(Matrix)
  library(MASS)
  library(lars,lib="/home/dun280/local/R-site/library")
  
  load(filename)
  #X<-scale(X,center=T,scale=F)
  #print(X[1:5,1:4])
  p<-dim(X)[2]
  nsamp<-dim(X)[1]
  
  foldslave <- function(folds){
    library(lars,lib="/home/dun280/local/R-site/library")
    library(Matrix)
    
    library(sparse.inv.cov,lib.loc="/home/dun280/local/R-site/library")
    
    get.neighbs<-function (X, kmax = trunc(nrow(X)/10), ff = 0, cols = 1:ncol(X), lars.type = "stepwise",
                           lars.normalise = TRUE, lars.intercept = TRUE, BIC01.only = FALSE,
                           prog_dir = "./", prog_filename = "get_neighbs.Rdata", prog_cycles = 100,
                           prog_exitonsave = FALSE, worker = "-1", restart = FALSE, restart_dir = "./")
      {
        require(Matrix)
        p <- ncol(X)
        N <- nrow(X)
        pen <- log(N)
        if (ff == 0)
          pen <- 2 * log(ncol(X))
        indx <- 1:p
        ii00 <- NULL
        ii05 <- NULL
        ii01 <- NULL
        ii2p <- NULL
        iiff <- NULL
        jj00 <- NULL
        jj05 <- NULL
        jj01 <- NULL
        jj2p <- NULL
        jjff <- NULL
        
        # The following if statement is to handle restarts after a checkpointed
        # job has exited early/erronously or finished their alotted number of cycles.
        # it first checks if a restart has been requested and whether the specified
        # restart directory contains a READABLE previous checkpoint file
        # for this worker  and dataset
        # if so, it overwrites the initial values of ii01 jj01 and cols appropriately
        # otherwise it lets processing start from the beginning
	if(restart){
          if(file.access(paste(restart_dir,"/",prog_filename,"_",worker,"_prog.Rdata", sep=""), mode=4) == 0)
            {
              load(paste(restart_dir,"/",prog_filename,"_",worker,"_prog.Rdata", sep=""))
              cols<-(i+1):cols[length(cols)]
              rm(i)
              gc()
            }
	}
        # write(paste("Worker",worker,"cols",cols[1],"-",cols[length(cols)],sep=" "),
        #  file=paste(prog_dir,"/temp_",worker,".Rdata", sep=""), append=TRUE)
        
	if(cols[1]<cols[length(cols)]){
          for (i in cols) {
            res <- reg.select(X[, -i], X[, i], kmax = kmax, ff = ff,
                              N = N, lars.type = lars.type, lars.normalise = lars.normalise,
                              lars.intercept = lars.intercept)
            gc()
            
            if (BIC01.only == TRUE) {
              vars <- (indx[-i])[res$selected[[3]]]
              q <- length(vars)
              if (q > 0) {
                ii01 <- c(ii01, rep(i, q))
                jj01 <- c(jj01, vars)
              }
              if( (i-cols[1]+1)%%prog_cycles == 0){
                # write(paste("Worker",worker,"saving at col",prog$cur,sep=" "),
                #file=paste(prog_dir,"/temp_",worker,".Rdata", sep=""), append=TRUE)
		save(ii01,jj01,i,file=paste(prog_dir,"/",prog_filename,"_",worker,"_prog.Rdata", sep=""))
		if(prog_exitonsave == TRUE){
                  return(list(a01="cp"))
		}
              }
              
            } # end (BIC01.only == TRUE)
            else {
              vars <- (indx[-i])[res$selected[[1]]]
              q <- length(vars)
              if (q > 0) {
                ii00 <- c(ii00, rep(i, q))
                jj00 <- c(jj00, vars)
              }
              vars <- (indx[-i])[res$selected[[2]]]
              q <- length(vars)
              if (q > 0) {
                ii05 <- c(ii05, rep(i, q))
                jj05 <- c(jj05, vars)
              }
              vars <- (indx[-i])[res$selected[[3]]]
              q <- length(vars)
              if (q > 0) {
                ii01 <- c(ii01, rep(i, q))
                jj01 <- c(jj01, vars)
              }
              vars <- (indx[-i])[res$selected[[4]]]
              q <- length(vars)
              if (q > 0) {
                ii2p <- c(ii2p, rep(i, q))
                jj2p <- c(jj2p, vars)
              }
              vars <- (indx[-i])[res$selected[[5]]]
              q <- length(vars)
              if (q > 0) {
                iiff <- c(iiff, rep(i, q))
                jjff <- c(jjff, vars)
              }
            }
          } # end  "for (i in cols)"
        } # end if(cols[1]<cols[length(cols)]){
        
        # Saving the final state is necessary so that the restart if-statment above
        # does not restart the calaultion by assuming no progress has been made
        #	prog<-list(ii01=ii01, jj01=jj01, i=cols[length(cols)])
	i<-cols[length(cols)]
	save(ii01, jj01, i,file=paste(prog_dir,"/",prog_filename,"_",worker,"_prog.Rdata", sep=""))
        
        
        if (BIC01.only == TRUE) {
          a01 <- new("dgTMatrix", i = as.integer(ii01 - 1),
                     j = as.integer(jj01 - 1), x = rep(1, length(ii01)), Dim = as.integer(c(p,p)))
          a01 <- as(a01, "dgCMatrix")
          a01 <- a01 + t(a01)
          diag(a01) <- 1
          return(list(a01 = a01))
        }
        else {
          a00 <- new("dgTMatrix", i = as.integer(ii00 - 1), j = as.integer(jj00 - 1),
                     x = rep(1, length(ii00)), Dim = as.integer(c(p,p)))
          
          a05 <- new("dgTMatrix", i = as.integer(ii05 - 1), j = as.integer(jj05 -1),
                     x = rep(1, length(ii05)), Dim = as.integer(c(p, p)))
          
          a01 <- new("dgTMatrix", i = as.integer(ii01 - 1), j = as.integer(jj01 - 1),
                     x = rep(1, length(ii01)), Dim = as.integer(c(p,p)))
          
          a2p <- new("dgTMatrix", i = as.integer(ii2p - 1), j = as.integer(jj2p - 1),
                     x = rep(1, length(ii2p)), Dim = as.integer(c(p, p)))
          
          aff <- new("dgTMatrix", i = as.integer(iiff - 1), j = as.integer(jjff - 1),
                     x = rep(1, length(iiff)), Dim = as.integer(c(p, p)))
          
          a00 <- as(a00, "dgCMatrix")
          a05 <- as(a05, "dgCMatrix")
          a01 <- as(a01, "dgCMatrix")
          a2p <- as(a2p, "dgCMatrix")
          aff <- as(aff, "dgCMatrix")
          a00 <- a00 + t(a00)
          a05 <- a05 + t(a05)
          a01 <- a01 + t(a01)
          a2p <- a2p + t(a2p)
          aff <- aff + t(aff)
          diag(a00) <- 1
          diag(a05) <- 1
          diag(a01) <- 1
          diag(a2p) <- 1
          diag(aff) <- 1
          return(list(a00 = a00, a05 = a05, a01 = a01, a2p = a2p,
                      aff = aff))
        }
      } #end get.neighbs
    
    
    #    Rprof(paste("./profile_",filename,"/slave_foldslave_", filename ,"_",mpi.comm.rank(),".profout", sep=""))
    
    load(filename)
    dd<-dim(X)
    #    paste(" dim(X) is ", dd[1], dd[2] , " on ", mpi.comm.rank(),"of",mpi.comm.size())
    
    #    paste(" 12. we are",mpi.comm.rank(),"of",mpi.comm.size())
    
    
    p<-dim(X)[2]
    remainder<- p%%folds
    
    if(foldNumber==folds){
      index <-1:floor(p/folds)+(foldNumber-1)*floor(p/folds)
      index<-c(index,c(max(index)+1):(max(index)+remainder))
      
    } else{
      index<-1:floor(p/folds)+(foldNumber-1)*floor(p/folds)
    }
    
    result<-get.neighbs(X,kmax = kmax,ff = 0, cols = index,
                        lars.type="stepwise",lars.normalise=FALSE, lars.intercept = TRUE,
                        BIC01.only=TRUE, prog_dir=outdir, prog_filename=filename,
                        prog_cycles=checkpointcycles, prog_exitonsave=exitoncheckpoint,
                        worker=mpi.comm.rank(), restart=restart, restart_dir=restartdir)
    
    #     result<-get.neighbs(X,kmax = kmax,ff = 0, cols = index,
    #                        lars.type="lar",lars.normalise=FALSE, lars.intercept = FALSE,
    #                      BIC01.only=TRUE)
    #    result<-get.neighbs(X,kmax = kmax,ff = 0, cols = index)
    
    #    Rprof(NULL)
    
    result$a01
    
    
  } # end foldslave
  
  mpi.bcast.Robj2slave(filename)
  mpi.bcast.Robj2slave(outdir)
  mpi.bcast.Robj2slave(restart)
  mpi.bcast.Robj2slave(restartdir)
  mpi.bcast.Robj2slave(checkpointcycles)
  mpi.bcast.Robj2slave(exitoncheckpoint)
  mpi.bcast.Robj2slave(kmax)
  
  mpi.bcast.Robj2slave(foldslave)
  
  
  result <- mpi.remote.exec(foldslave(mpi.comm.size()-1))
  
  save_file<-TRUE
  count<-0
  for ( i in result ){
    if(class(i)== "character"){
      if( i=="cp"){
        save_file<-FALSE
        print("one of the workers returned checkpoint")
      }
      print("worker returned a character string, but not cp")
    }
    count<-count+1
    print(paste("worker ",count,"s results checked by master",sep=""))
  }
  
  if(save_file){
    print("saving result file!")
    #  	temp3 <- Reduce("+", result) #commented out. I think the Reduce may be having problems with
    #  	temp3[temp3!=0]<-1           # the size of the matrix
    #  
    #  	if (mpi.comm.rank() == 0) {
    #    		save(temp3,file=paste("result",filename, sep="_"))
    #  	}
  } else {
    print("At least one slave checkpointed but exited with out finishing")
    print(paste("Job should be restarted using the chechpoint files in",outdir,sep=" "))
  }
  
  mpi.close.Rslaves()
  Sys.sleep(5)
  mpi.quit()
  
}

#Rprof(paste("./profile_",filename,"/master_myfunc_", filename, ".profout", sep=""))
myfunc()
#Rprof(NULL)

