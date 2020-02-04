#' Berger-Parker computation function
#' 
#' Compute the Berger-Parker index, after data check
#' 
#' @param input the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it
#' @param window the size of the square window size. Default value is 3.
#' @param simplify the 10 power which will be used to covert float into natural number. Default value is 3.
#' @param nc.cores the nuber of cores which will be used. Default value is 1.
#' @param cluster.type the type of cluster which will be used. Default type is "MPI".
#' @param debugging a boolean variable set to FALSE by default. If TRUE, let the user check all the steps
#' 
#' @return Matrix or a list of matrixes with the Berger-Parker index computed through moving window of the given size
#'
#'
BergerParker <- function(input, window=3, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
    # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
    if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100 ^ simplify
      rasterm <- getValues(rasterm) * mfactor
      rasterm <- as.integer(rasterm)
      rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
      gc()
    }
    else{
      rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
    }
  }
  #Print user messages
  else if( is(input,"matrix") | is(input,"list") ) {
    isfloat<-FALSE # If data are float numbers, transform them in integer
    if( !is.integer(rasterm) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100^simplify
      rasterm <- as.integer(rasterm*mfactor)
      rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
      gc()
    }else{
      rasterm<-as.matrix(rasterm)
    }
  }
  message("Matrix check OK: \nBerger-Parker output matrix will be returned")
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    outS <- BergerParkerS(rasterm, w, debugging)
    message(("\nCalculation of Berger-Parker's index is complete!\n"))
    return (outS)
  }
  else {
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Berger-Parker parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    outP <- BergerParkerP(rasterm, w,  debugging)
    outP <- do.call(cbind,outP)
    return(outP)
  }
}
#' Berger-Parker single core computation function
#' 
#' Compute the Berger-Parker index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A matrix with the Berger-Parker index computed through moving window associated to
#'
#'


BergerParkerS <- function(rasterm, w,  debugging){
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  message("\nStarting Berger-Parker index calculation:\n")
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add "fake" columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Berger-Parker\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_values<-as.vector(tw)
      p<-max(tw_values/sum(tw_values))
      out[rw-w,cl-w]<-(log(1/p))
    }   
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } 
  return(out)
}
#' Berger-Parker multiple cores computation function
#' 
#' Compute the Berger-Parker index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A list of vectors with the Berger-Parker index computed through moving window associated to
#'
#'

BergerParkerP<-function(rasterm, w, debugging){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  BergerParkerOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    BergerParkerOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Berger-Parker - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- max(tw_values/sum(tw_values))
      vv <- (log(1/p))
    })
    return(BergerParkerOut)
  } # End Berger-Parker - parallelized
  message(("\n\n Parallel calculation of Berger-Parker's index complete.\n"))
  return(BergerParkerOP)
}

#' Shannon computation function
#' 
#' Compute the Shannon index, after data check
#' 
#' @param input the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.
#' @param window the size of the square window size. Default value is 3.
#' @param simplify the 10 power which will be used to covert float into natural number. Default value is 3.
#' @param nc.cores the nuber of cores which will be used. Default value is 1.
#' @param cluster.type the type of cluster which will be used. Default type is "MPI".
#' @param debugging a boolean variable set to FALSE by default. If TRUE, let the user check all the steps
#' 
#' @return Matrix or a list of matrixes with the Shannon index computed through moving window of the given size
#'
#'
Shannon <- function(input, type, window=3, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
    if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100 ^ simplify
      rasterm <- getValues(rasterm) * mfactor
      rasterm <- as.integer(rasterm)
      rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
      gc()
    }
    else{
      rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
    }
  }
  else if( is(input,"matrix") | is(input,"list") ) {
    isfloat<-FALSE # If data are float numbers, transform them in integer
    if( !is.integer(rasterm) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100^simplify
      rasterm <- as.integer(rasterm*mfactor)
      rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
      gc()
    }else{
      rasterm<-as.matrix(rasterm)
    }
  }
  #Print user messages
  message("Matrix check OK: \nShannon output matrix will be returned")
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    outS <- ShannonS(rasterm, w, debugging)
    message(("\nCalculation of Shannon's index is complete!\n"))
    return (outS)
  }
  if (nc.cores>1){
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Renyi parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    outP <- ShannonP(rasterm, w, debugging)
    outP <- do.call(cbind,outP)
    return(outP)
  }
}

#' Shannon single core computation function
#' 
#' Compute the Shannon index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A matrix with the Shannon index computed through moving window associated to
#'
#'

ShannonS <- function(rasterm, w, debugging){
  message("\nStarting Shannon-Wiener index calculation:\n")
  # Reshape values
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add "fake" columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Shannon-Wiener\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_values<-as.vector(tw)
      p<-tw_values/sum(tw_values)
      p_log<-log(p)
      out[rw-w,cl-w]<-(-(sum(p*p_log)))
      #}
    }   
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } 
  
  return(out)
}

#' Shannon multiple cores computation function
#' 
#' Compute the Shannon index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A list of vectors with the Shannon index computed through moving window associated to
#'
#'

ShannonP<-function(rasterm, w, debugging){
  #
  ## Reshape values
  #
  values <- as.numeric( as.factor(rasterm) )
  rasterm_1 <- matrix(data = values, nrow = dim(rasterm)[1], ncol = dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = dim(rasterm)[2], nrow = w)
  ver <- matrix(NA, ncol = w, nrow = dim(rasterm)[1]+ w * 2)
  trasterm <- cbind(ver, rbind(hor,rasterm_1,hor), ver)
  rm(hor, ver, rasterm_1, values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ShannonOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    ShannonOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Shannon- parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      vv <- (-(sum(p*log(p))))
      return(vv)
      #}
    })
    return(ShannonOut)
  } # End Shannon- parallelized
  message(("\n\n Parallel calculation of Shannon's index complete.\n"))
  return(ShannonOP)
}

#' Pielou computation function
#' 
#' Compute the Pielou index, after data check
#' 
#' @param input the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.
#' @param window the size of the square window size. Default value is 3.
#' @param simplify the 10 power which will be used to covert float into natural number. Default value is 3.
#' @param nc.cores the nuber of cores which will be used. Default value is 1.
#' @param cluster.type the type of cluster which will be used. Default type is "MPI".
#' @param debugging a boolean variable set to FALSE by default. If TRUE, let the user check all the steps
#' 
#' @return Matrix or a list of matrixes with the Pielou index computed through moving window of the given size
#'
#'

Pielou <- function(input, window=3, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
    if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100 ^ simplify
      rasterm <- getValues(rasterm) * mfactor
      rasterm <- as.integer(rasterm)
      rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
      gc()
    }
    else{
      rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
    }
  }
  else if( is(input,"matrix") | is(input,"list") ) {
    isfloat<-FALSE # If data are float numbers, transform them in integer
    if( !is.integer(rasterm) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100^simplify
      rasterm <- as.integer(rasterm*mfactor)
      rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
      gc()
    }else{
      rasterm<-as.matrix(rasterm)
    }
  }
  #Print user messages
  message("Matrix check OK: \nPielou output matrix will be returned")
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    outS <- PielouS(rasterm, w, debugging)
    message(("\nCalculation of Pielou's index is complete!\n"))
    return (outS)
  }
  if (nc.cores>1){
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Pielou parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    outP <- PielouP(rasterm, w, debugging)
    return(outP)
  }
}

#' Pielou single core computation function
#' 
#' Compute the Pielou index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A matrix with the Pielou index computed through moving window associated to
#'
#'


PielouS <- function(rasterm, w, debugging){
  message("\nStarting Pielou's index calculation:\n")
  # Reshape values
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add "fake" columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  window <- 2*w + 1
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Pielou\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_values<-as.vector(tw)
      p<-tw_values/sum(tw_values)
      p_log<-log(p)
      out[rw-w,cl-w]<-(-(sum(p*p_log))/log((window^2)))
      #}
    }   
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } 
  
  return(out)
}


#' Pielou multiple cores computation function
#' 
#' Compute the Pielou index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A list of vectors with the Pielou index computed through moving window associated to
#'
#'

PielouP<-function(rasterm, w, debugging){
  #
  ## Reshape values
  #
  values <- as.numeric( as.factor(rasterm) )
  rasterm_1 <- matrix(data = values, nrow = dim(rasterm)[1], ncol = dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = dim(rasterm)[2], nrow = w)
  ver <- matrix(NA, ncol = w, nrow = dim(rasterm)[1]+ w * 2)
  trasterm <- cbind(ver, rbind(hor,rasterm_1,hor), ver)
  rm(hor, ver, rasterm_1, values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  window <- 2*w + 1
  PielouOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    PielouOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Pielou - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      vv <- (-(sum(p*log(p)))/log((window^2)))
      return(vv)
      #}
    })
    return(PielouOut)
  } # End Pielou - parallelized
  message(("\n\n Parallel calculation of Pielou's index complete.\n"))
  return(PielouOP)
}

#' Rényi computation function
#' 
#' Compute the Rényi index, after data check
#' 
#' @param input the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.
#' @param window the size of the square window size. Default value is 3.
#' @param mode a variable which can be 
#' \itemize{
#' \item "single" to compute Rényi index for just one alpha value, the dealt value
#' \item "iterative" to compute Rényi index for all the integers value of alpha in a given interval,
#' \item "sequential" to compute Rényi index for all the alpha values in a given vector.
#' }
#' @param BergerParker a boolean variable. The default value is FALSE. If TRUE, in the "single" mode the function will return just the Berger-Parker index. If TRUE, in all the other modes, it will return also the Berger-Parker index.
#' @param simplify the 10 power which will be used to covert float into natural number. Default value is 3.
#' @param alpha a numerical variable. Its default value is 1.
#' \itemize{
#' \item In "single" mode, alpha has to be a numerical value greater than 0;
#' \item In "iterative" mode, alpha has to be a length 2 vector;
#' \item In "sequential" mode, alpha has to be a vector of leangth at least 2.
#' }
#' @param base a numerical value, which lets the user choose the base of the logarithm in Rényi index formula. Its defalut value is exp(1) 
#' @param nc.cores the nuber of cores which will be used. Default value is 1.
#' @param cluster.type the type of cluster which will be used. Default type is "MPI".
#' @param debugging a boolean variable set to FALSE by default. If TRUE, let the user check all the steps
#' 
#' @return Matrix or a list of matrixes with the Rényi index computed through moving window of the given size
#'
#'

Renyi <- function(input, window=3, mode ="single", BergerParker=FALSE, alpha=1, base=exp(1), simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  Shannon <- FALSE
  Hill <- FALSE
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  if (base<.Machine$double.eps){
    stop("base value must be in the (0,+\u221E) interval. Exiting...")
  }
  if (!is.numeric(alpha)){
    stop("alpha must be a number or a vector. Exiting...")
  }
  if (mode=="single"){
    if (length(alpha)!=1){
      stop("In mode \"single\" alpha must be a single number. Exiting...")
    }
    if (alpha<0){
      stop("The alpha value must be a non-negative number. Exiting...")
    }
    if(abs(alpha-1)<.Machine$double.eps && !BergerParker){
      Shannon <- TRUE
    }
    if(alpha >=.Machine$integer.max){
      BergerParker <- TRUE
    }
    if(integer){
      alpha <- as.integer(alpha)
    }
  }
  else if(mode == "iterative"){
    if (length(alpha) != 2){
      stop("In mode \"iterative\" alpha must be a numeric vector containing the star and the stop value. Exiting...")
    }
    start <- as.integer(alpha[1])
    if (start<0){
      stop("The starting value must be a non-negative number. Exiting...")
    }
    stop <- as.integer(alpha[2])
    if (stop <= start){
      stop("Integer part of the starting value, alpha[1], must be strictly greater that the integer part of the stopping value, alpha[2]. Exiting...")
    }
    if (start <= 1 && 1 <= stop ){
      Shannon <- TRUE
    }
  }
  else if(mode == "sequential"){
    if ( length(alpha) < 2){
      stop("In mode \"sequential\" alpha must be a numeric vector containing at least two values. Exiting...")
    }
    if ( length(which(alpha < 0)) != 0 ){
      stop("The alpha values must be non-negative numbers. Exiting...")
    }
    if ( integer ){
      val <- unique( as.integer(alpha) )
    }
    else {
      val <- unique(alpha)
    }
  }
  else{
    stop("The choosen mode is not defined. Exiting...")
  }
  
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
    if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100 ^ simplify
      rasterm <- getValues(rasterm) * mfactor
      rasterm <- as.integer(rasterm)
      rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
      gc()
    }
    else{
      rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
    }
    #Print user messages
    if (mode == "single"){
      if( BergerParker ){
        message("Matrix check OK: \nBerger-Parker output matrix will be returned")
      }else if( Shannon ){
        message("Matrix check OK: \nShannon output matrix will be returned")
      }else {
        message("Matrix check OK: \nRenyi with parameter value=", alpha," output matrix will be returned")
      }
    }
    else if (mode == "iterative"){
      if( BergerParker && !Shannon){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
      }else if( BergerParker && Shannon  ){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
      }else if( !BergerParker && Shannon  ){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
      }else {
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"]")
      }
    }
    # If data are a matrix or a list
  }else if( is(input,"matrix") | is(input,"list") ) {
    isfloat<-FALSE # If data are float numbers, transform them in integer
    if( !is.integer(rasterm) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100^simplify
      rasterm <- as.integer(rasterm*mfactor)
      rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
      gc()
    }else{
      rasterm<-as.matrix(rasterm)
    }
    #Print user messages
    if (mode == "single"){
      if( BergerParker ){
        message("Matrix check OK: \nBerger-Parker output matrix will be returned")
      }else if( Shannon ){
        message("Matrix check OK: \nShannon output matrix will be returned")
      }else {
        message("Matrix check OK: \nRenyi with parameter value=", alpha," output matrix will be returned")
      }
    }
    else if (mode == "iterative"){
      if( BergerParker && !Shannon){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
      }else if( BergerParker && Shannon  ){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
      }else if( !BergerParker && Shannon  ){
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
      }else {
        message("Matrix check OK: \nRenyi output matrix will be returned for parameters integer values in [",start,",",stop,"]")
      }
    }
  }
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    if(mode == "single") {
      if( Shannon ) {
        outS <- ShannonS(rasterm, w, debugging)
        message(("\nCalculation of Shannon's index is complete!\n"))
      } # End ShannonD
      else if( BergerParker ) {
        outS <- BergerParkerS(rasterm, w, debugging)
        message(("\nCalculation of Berger-Parker's index is complete!\n"))
      } # End BergerParker
      else{
        outS<- RenyiS(rasterm, w, alpha, base, debugging)
        message(("\nCalculation of Renyi's index complete.\n"))
      }
      return (outS)
    }
    else if(mode == "iterative"){
      out <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out[[s]] <- ShannonS(rasterm, w, debugging)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RenyiAlpha",as.character(ALPHA))
          out[[s]] <- RenyiS(rasterm,w, ALPHA, base, debugging)
        }
      }
      message(("\nCalculation of Renyi's index complete.\n"))
      if (BergerParker){
        s<-"Berger-Parker"
        out[[s]] <- BergerParkerS(rasterm, w, debugging)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      
      return(out)
    }
    else if(mode == "sequential"){
      out <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out[[s]] <- ShannonS(rasterm, w, debugging)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RenyiAlpha",as.character(ALPHA))
          out[[s]] <- RenyiS(rasterm, w, ALPHA, base, debugging)
        }
      }
      message(("\nCalculation of Renyi's index complete.\n"))
      if (BergerParker){
        s<-"Berger-Parker"
        out[[s]] <- BergerParkerS(rasterm, w, debugging)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      
      return(out)
    }
    
  }
  if (nc.cores>1){
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Renyi parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    
    if(mode == "single") {
      if( Shannon ) {
        outP <- ShannonP(rasterm, w, debugging)
      }
      else if (BergerParker){
        outP <- BergerParkerP(rasterm, w, debugging)
      }
      else{
        outP <- RenyiP(rasterm, w, alpha, base, debugging)
      }
      return(do.call(cbind,outP))
    }
    else if(mode == "iterative"){
      outP <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, debugging)
          outP[[s]] <- do.call(cbind,out)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RenyiAlpha",as.character(ALPHA))
          out<- RenyiP(rasterm, w, ALPHA, base, debugging)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Renyi's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, debugging)
        outP[[s]] <- do.call(cbind,out)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      return(outP)
    }
    else if( mode == "sequential"){
      outP <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, debugging)
          outP[[s]] <- do.call(cbind,out)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RenyiAlpha",as.character(ALPHA))
          out<- RenyiP(rasterm, w, ALPHA, base, debugging)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Renyi's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, debugging)
        outP[[s]] <- do.call(cbind,out)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      return(outP)
    }
  }
}

#' Rényi single core computation function
#' 
#' Compute the Rényi index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param alpha a numerical variable
#' @param base a numerical value, which lets the user choose the base of the logarithm in Rényi index formula. Its defalut value is exp(1) 
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A matrix with the Rényi index computed with the fixed alpha through moving window associated to w.
#'

RenyiS <- function(rasterm, w, alpha, base, debugging){
  message("\n\nStarting Renyi index calculation with parameter value = ",alpha," \n\n\n\n\n")
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add fake columns and rows for moving window
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Renyi\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      #if clause to exclude windows with only 1 category
      #if(length(tw_values) < 2) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      p <- tw_values/sum(tw_values)
      out[rw-w,cl-w]<-1/(1-alpha) * drop(log(sum(p^alpha),base))
      #}
      #} 
    }
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } # End of for loop 
  return(out)
}

#' Rényi multiple cores computation function
#' 
#' Compute the Rényi index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param alpha a numerical variable
#' @param base a numerical value, which lets the user choose the base of the logarithm in Rényi index formula. Its defalut value is exp(1) 
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A list of vectors with the Rényi index computed with the fixed alpha through moving window associated to w.
#'


RenyiP<-function(rasterm, w, alpha, base, debugging){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  RenyiOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    RenyiOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Renyi - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      vv <- 1/(1-alpha) * drop(log(sum(p^alpha),base))
      return(vv)
      #}
    })
    return(RenyiOut)
  } # End classic Renyi - parallelized
  message(("\n\n Parallel calculation of Renyi's index complete.\n"))
  return(RenyiOP)
}


#' Hill computation function
#' 
#' Compute the Hill index, after data check
#' 
#' @param input the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.
#' @param window the size of the square window size. Default value is 3.
#' @param mode a variable which can be 
#' \itemize{
#' \item "single" to compute Hill index for just one alpha value, the dealt value
#' \item "iterative" to compute Hill index for all the integers value of alpha in a given interval,
#' \item "sequential" to compute Hill index for all the alpha values in a given vector.
#' }
#' @param BergerParker a boolean variable. The default value is FALSE. If TRUE, in the "single" mode the function will return just the exponential of Berger-Parker index. If TRUE, in all the other modes, it will return also the exponential of Berger-Parker Berger-Parker index.
#' @param simplify the 10 power which will be used to covert float into natural number. Default value is 3.
#' @param alpha a numerical variable. Its default value is 1.
#' \itemize{
#' \item In "single" mode, alpha has to be a numerical value greater than 0;
#' \item In "iterative" mode, alpha has to be a length 2 vector;
#' \item In "sequential" mode, alpha has to be a vector of leangth at least 2.
#' }
#' @param nc.cores the nuber of cores which will be used. Default value is 1.
#' @param cluster.type the type of cluster which will be used. Default type is "MPI".
#' @param debugging a boolean variable set to FALSE by default. If TRUE, let the user check all the steps
#' 
#' @return Matrix or a list of matrixes with the Hill index computed through moving window of the given size
#'
#'


Hill <- function(input, window=3, mode ="single", BergerParker=FALSE, alpha=1, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  Shannon <- FALSE
  Hill <- FALSE
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  if (!is.numeric(alpha)){
    stop("alpha must be a number or a vector. Exiting...")
  }
  if (mode=="single"){
    if (length(alpha)!=1){
      stop("In mode \"single\" alpha must be a single number. Exiting...")
    }
    if (alpha<0){
      stop("The alpha value must be a non-negative number. Exiting...")
    }
    if(abs(alpha-1)<.Machine$double.eps && !BergerParker){
      Shannon <- TRUE
    }
    if(alpha >=.Machine$integer.max){
      BergerParker <- TRUE
    }
    if(integer){
      alpha <- as.integer(alpha)
    }
  }
  else if(mode == "iterative"){
    if (length(alpha) != 2){
      stop("In mode \"iterative\" alpha must be a numeric vector containing the star and the stop value. Exiting...")
    }
    start <- as.integer(alpha[1])
    if (start<0){
      stop("The starting value must be a non-negative number. Exiting...")
    }
    stop <- as.integer(alpha[2])
    if (stop <= start){
      stop("Integer part of the starting value, alpha[1], must be strictly greater that the integer part of the stopping value, alpha[2]. Exiting...")
    }
    if (start <= 1 && 1 <= stop ){
      Shannon <- TRUE
    }
  }
  else if(mode == "sequential"){
    if ( length(alpha) < 2){
      stop("In mode \"sequential\" alpha must be a numeric vector containing at least two values. Exiting...")
    }
    if ( length(which(alpha < 0)) != 0 ){
      stop("The alpha values must be non-negative numbers. Exiting...")
    }
    if ( integer ){
      val <- unique( as.integer(alpha) )
    }
    else {
      val <- unique(alpha)
    }
  }
  else{
    stop("The choosen mode is not defined. Exiting...")
  }
  
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
    if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100 ^ simplify
      rasterm <- getValues(rasterm) * mfactor
      rasterm <- as.integer(rasterm)
      rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
      gc()
    }
    else{
      rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
    }
    #Print user messages
    if (mode == "single"){
      if( BergerParker ){
        message("Matrix check OK: \nBerger-Parker output matrix will be returned")
      }else if( Shannon ){
        message("Matrix check OK: \nShannon output matrix will be returned")
      }else {
        message("Matrix check OK: \nHill with parameter value=", alpha," output matrix will be returned")
      }
    }
    else if (mode == "iterative"){
      if( BergerParker && !Shannon){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
      }else if( BergerParker && Shannon  ){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
      }else if( !BergerParker && Shannon  ){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
      }else {
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"]")
      }
    }
    # If data are a matrix or a list
  }else if( is(input,"matrix") | is(input,"list") ) {
    isfloat<-FALSE # If data are float numbers, transform them in integer
    if( !is.integer(rasterm) ){
      message("Converting input data in an integer matrix...")
      isfloat <- TRUE
      mfactor <- 100^simplify
      rasterm <- as.integer(rasterm*mfactor)
      rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
      gc()
    }else{
      rasterm<-as.matrix(rasterm)
    }
    #Print user messages
    if (mode == "single"){
      if( BergerParker ){
        message("Matrix check OK: \nBerger-Parker output matrix will be returned")
      }else if( Shannon ){
        message("Matrix check OK: \nShannon output matrix will be returned")
      }else {
        message("Matrix check OK: \nHill with parameter value=", alpha," output matrix will be returned")
      }
    }
    else if (mode == "iterative"){
      if( BergerParker && !Shannon){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
      }else if( BergerParker && Shannon  ){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
      }else if( !BergerParker && Shannon  ){
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
      }else {
        message("Matrix check OK: \nHill output matrix will be returned for parameters integer values in [",start,",",stop,"]")
      }
    }
  }
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    if(mode == "single") {
      if( Shannon ) {
        outS <- exp(ShannonS(rasterm, w, debugging))
        message(("\nCalculation of Shannon's index is complete!\n"))
      } # End ShannonD
      else if( BergerParker ) {
        outS <- exp(BergerParkerS(rasterm, w, debugging))
        message(("\nCalculation of Berger-Parker's index is complete!\n"))
      } # End BergerParker
      else{
        outS<- HillS(rasterm, w, alpha, debugging)
        message(("\nCalculation of Hill's index complete.\n"))
      }
      return (outS)
    }
    else if(mode == "iterative"){
      out <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out[[s]] <- exp(ShannonS(rasterm, w, debugging))
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("HillAlpha",as.character(ALPHA))
          out[[s]] <- HillS(rasterm,w, ALPHA, debugging)
        }
      }
      message(("\nCalculation of Hill's index complete.\n"))
      if (BergerParker){
        s<-"Berger-Parker"
        out[[s]] <- exp(BergerParkerS(rasterm, w, debugging))
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      
      return(out)
    }
    else if(mode == "sequential"){
      out <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out[[s]] <- exp(ShannonS(rasterm, w, debugging))
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("HillAlpha",as.character(ALPHA))
          out[[s]] <- HillS(rasterm, w, ALPHA, debugging)
        }
      }
      message(("\nCalculation of Hill's index complete.\n"))
      if (BergerParker){
        s<-"Berger-Parker"
        out[[s]] <- exp(BergerParkerS(rasterm, w, debugging))
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      
      return(out)
    }
    
  }
  if (nc.cores>1){
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Hill parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    
    if(mode == "single") {
      if( Shannon ) {
        outP <- exp(ShannonP(rasterm, w, debugging))
      }
      else if (BergerParker){
        outP <- exp(BergerParkerP(rasterm, w, debugging))
      }
      else{
        outP <- HillP(rasterm, w, alpha, debugging)
      }
      return(do.call(cbind,outP))
    }
    else if(mode == "iterative"){
      outP <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, debugging)
          outP[[s]] <- exp(do.call(cbind,out))
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("HillAlpha",as.character(ALPHA))
          out<- HillP(rasterm, w, ALPHA, debugging)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Hill's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, debugging)
        outP[[s]] <- exp(do.call(cbind,out))
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      return(outP)
    }
    else if( mode == "sequential"){
      outP <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, debugging)
          outP[[s]] <- exp(do.call(cbind,out))
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("HillAlpha",as.character(ALPHA))
          out<- HillP(rasterm, w, ALPHA, debugging)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Hill's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, debugging)
        outP[[s]] <- exp(do.call(cbind,out))
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      return(outP)
    }
  }
}

#' Hill single core computation function
#' 
#' Compute the Hill index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param alpha a numerical variable
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A matrix with the Hill index computed with the fixed alpha through moving window associated to w.
#'

HillS <- function(rasterm, w, alpha, debugging){
  message("\n\nStarting Hill index calculation with parameter value = ",alpha," \n\n\n\n\n")
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add fake columns and rows for moving window
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Hill\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      #if clause to exclude windows with only 1 category
      #if(length(tw_values) < 2) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      p <- tw_values/sum(tw_values)
      out[rw-w,cl-w] <- drop(sum(p^alpha))^(1/(1-alpha)) 
      #}
      #} 
    }
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } # End of for loop 
  return(out)
}

#' Hill multiple cores computation function
#' 
#' Compute the Hill index, without any check over the data
#' 
#' @param rasterm a matrix, over which the index will be computed.
#' @param w a value obatined from the window side, subtracting 1 and dividing by 2.
#' @param alpha a numerical variable
#' @param debugging a boolean variable. If TRUE, let the user check all the steps.
#' 
#' @return A list of vectors with the Hill index computed with the fixed alpha through moving window associated to w.
#'

HillP<-function(rasterm, w, alpha, debugging){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  HillOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    HillOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Hill - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      vv <- drop(sum(p^alpha))^(1/(1-alpha))
      return(vv)
      #}
    })
    return(HillOut)
  } # End classic Hill - parallelized
  message(("\n\n Parallel calculation of Hill's index complete.\n"))
  return(HillOP)
}


######### SPECTRALRAO #############################
## Developed by Matteo Marcantonio
## Latest update: 18th March 2019
## -------------------------------------------------
## Code to calculate Rao's quadratic entropy on a
## numeric matrix, RasterLayer object (or lists)
## using a moving window algorithm. 
## The function also calculates Shannon-Wiener index.
## -------------------------------------------------
## Rao's Q Min = 0, if all pixel classes have
## distance 0. If the chosen distance ranges between
## 0 and 1, Rao's Max = 1-1/S (Simpson Diversity,
## where S is the number of pixel classes).
## -------------------------------------------------
## Find more info and application here: 
## 1) https://doi.org/10.1016/j.ecolind.2016.07.039 
## 2) https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12941
####################################################
# Function
spectralrao <- function(input, distance_m="euclidean", p=NULL, window=9, mode="classic", lambda=0, shannon=FALSE, rescale=FALSE, na.tolerance=0.0, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE, ...)
{
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  require(proxy)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  if(na.tolerance>1){
    stop("na.tolerance must be in the [0-1] interval. Exiting...")
  }
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
    if( mode=="classic" ){
      isfloat<-FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
      if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
        message("Converting input data in an integer matrix...")
        isfloat<-TRUE
        mfactor<-100^simplify
        rasterm<-getValues(rasterm)*mfactor
        rasterm<-as.integer(rasterm)
        rasterm<-matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
        gc()
      }else{
        rasterm<-matrix(getValues(rasterm),ncol=ncol(input),nrow=nrow(input),byrow=TRUE)
      }
    }
    #Print user messages
    if( mode=="classic" & shannon ){
      message("Matrix check OK: \nRao and Shannon output matrices will be returned")
    }else if( mode=="classic" & !shannon ){
      message("Matrix check OK: \nRao output matrix will be returned")
    }else if( mode=="multidimension" & !shannon ){
      message(("Matrix check OK: \nA matrix with multimension RaoQ will be returned"))
    }else if( mode=="multidimension" & shannon ){
      stop("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE")
    }else{
      stop("Matrix check failed: \nNot a valid input | method | distance, please check all these options...")
    }
    # If data are a matrix or a list
  }else if( is(input,"matrix") | is(input,"list") ) {
    if( mode=="classic" ){ 
      isfloat<-FALSE # If data are float numbers, transform them in integer
      if( !is.integer(rasterm) ){
        message("Converting input data in an integer matrix...")
        isfloat<-TRUE
        mfactor<-100^simplify
        rasterm<-as.integer(rasterm*mfactor)
        rasterm<-matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
        gc()
      }else{
        rasterm<-as.matrix(rasterm)
      }
    }
    if( mode=="classic" & shannon ){
      message("Matrix check OK: \nRao and Shannon output matrices will be returned")
    }else if( mode=="classic" & !shannon ){
      message("Matrix check OK: \nRao output matrix will be returned")
    }else if( mode=="multidimension" & shannon ){
      stop("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE")
    }else if( mode=="multidimension" & !shannon ){
      message(("Matrix check OK: \nA matrix with multimension RaoQ will be returned"))
    }else{
      stop("Matrix check failed: \nNot a valid input | method | distance, please check all these options")
    }
  }
  
  if(nc.cores>1) {
    if(mode=="multidimension"){
      message(
        "Multi-core is not supported for multidimensional Rao, proceeding with 1 core...")
      nc.cores=1
    }else{
      message("
##################### Starting parallel calculation #######################")
    }
  }
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  #
  ## Preparation of output matrices
  #
  if(nc.cores==1) {
    raoqe<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  }
  if(shannon){
    shannond<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  }
  #
  ## If mode is classic Rao
  #
  if(mode=="classic") {
    #
    # If classic RaoQ is parallelized
    #
    if(nc.cores>1) {
      #
      ## Required packages for parallel calculation
      #
      require(foreach)
      require(doSNOW)
      require(parallel)
      if( cluster.type=="MPI" ){
        require(Rmpi)
      }
      #
      ## Reshape values
      #
      values<-as.numeric(as.factor(rasterm))
      rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
      #
      ## Add additional columns and rows to match moving window
      #
      hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
      ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
      trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
      rm(hor,ver,rasterm_1,values); gc()
      if(debugging){cat("#check: RaoQ parallel function.")}
      #       
      ## Derive distance matrix
      #
      if( is.character( distance_m) | is.function(distance_m) ) {
        d1<-proxy::dist(as.numeric(levels(as.factor(rasterm))),method=distance_m)
      } else if( is.matrix(distance_m) | is.data.frame(distance_m) ) {
        d1<-stats::as.dist(xtabs(distance_m[, 3] ~ distance_m[, 2] + distance_m[, 1]))
      }
      #
      ## Export variables in the global environment
      #
      if(isfloat) {
        sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
      }
      #       
      ## Create cluster object with given number of slaves
      #
      plr<<-TRUE
      if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
        cls <- parallel::makeCluster(nc.cores,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
      } else if( cluster.type=="MPI" ) {
        cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
      }
      registerDoSNOW(cls)
      clusterCall(cl=cls, function() library("parallel"))
      if(isfloat) {
        parallel::clusterExport(cl=cls, varlist=c("mfactor"))
      }
      on.exit(stopCluster(cls)) # Close the clusters on exit
      gc()
      #
      ## Start the parallelized loop over iter
      #
      pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      raop <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
        if(debugging) {
          cat(paste(cl))
        }
        raout <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
          if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
            vv<-NA
            return(vv)
          } 
          else {
            tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
              tw<-tw[-length(tw)]
            }
            if( debugging ) {
              message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
            if( length(tw_values) <2 ) {
              vv<-NA
              return(vv)
            }
            else {
              p <- tw_values/sum(tw_values)
              p1 <- diag(0,length(tw_values))
              p1[upper.tri(p1)] <- c(combn(p,m=2,FUN=prod))
              d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
              vv <- sum(p1*d2)
              return(vv)
            }
          }
        })
        return(raout)
      } # End classic RaoQ - parallelized
      message(("\n\nCalculation of Rao's index complete.\n"))
      #
      ## If classic RaoQ is sequential
      #
    } else if(nc.cores==1) {
      # Reshape values
      values<-as.numeric(as.factor(rasterm))
      rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
      # Add fake columns and rows for moving window
      hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
      ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
      trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
      # Derive distance matrix
      classes<-levels(as.factor(rasterm))
      if( is.character(distance_m) | is.function(distance_m) ) {
        d1<-proxy::dist(as.numeric(classes),method=distance_m)
      } else if( is.matrix(distance_m) | is.data.frame(distance_m) ) {
        d1<-stats::as.dist(xtabs(distance_m[, 3] ~ distance_m[, 2] + distance_m[, 1]))
      }
      # Loop over each pixel
      for (cl in (1+w):(dim(rasterm)[2]+w)) {
        for(rw in (1+w):(dim(rasterm)[1]+w)) {
          if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
            raoqe[rw-w,cl-w]<-NA
          } else {
            tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
              tw<-tw[-length(tw)]
            }
            if(debugging) {
              message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
            if(length(tw_values) < 2) {
              raoqe[rw-w,cl-w]<-NA
            } else {
              p <- tw_values/sum(tw_values)
              p1 <- diag(0,length(tw_values))
              p1[upper.tri(p1)] <- c(combn(p,m=2,FUN=prod))
              d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
              if(isfloat) {
                raoqe[rw-w,cl-w]<-sum(p1*d2)/mfactor
              } else {
                raoqe[rw-w,cl-w]<-sum(p1*d2)
              }
            }
          } 
          progress(value=cl, max.value=c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2, progress.bar = FALSE)
        } 
      } # End of for loop 
      message(("\nCalculation of Rao's index complete.\n"))
    }
  }  # End classic RaoQ - sequential
  else if( mode=="multidimension" ){
    if(debugging) {
      message("#check: Into multidimensional clause.")
    }
    #----------------------------------------------------#
    #
    ## If multimensional RaoQ
    #
    # Check if there are NAs in the matrices
    if ( is(rasterm,"RasterLayer") ){
      if(any(sapply(lapply(unlist(input),length),is.na)==TRUE))
        message("\n Warning: One or more RasterLayers contain NA which will be threated as 0")
    } else if ( is(rasterm,"matrix") ){
      if(any(sapply(input, is.na)==TRUE) ) {
        message("\n Warning: One or more matrices contain NA which will be threated as 0")
      }
    }
    #
    ## Check whether the chosen distance metric is valid or not
    #
    if( distance_m=="euclidean" | distance_m=="manhattan" | distance_m=="canberra" | distance_m=="minkowski" | distance_m=="mahalanobis" ) {
      #
      ## Define distance functions
      #
      #euclidean
      multieuclidean <- function(x) {
        tmp <- lapply(x, function(y) {
          (y[[1]]-y[[2]])^2
        })
        return(sqrt(Reduce(`+`,tmp)))
      }
      #manhattan
      multimanhattan <- function(x) {
        tmp <- lapply(x, function(y) {
          abs(y[[1]]-y[[2]])
        })
        return(Reduce(`+`,tmp))
      }
      #canberra
      multicanberra <- function(x) {
        tmp <- lapply(x, function(y) {
          abs(y[[1]] - y[[2]]) / (abs(y[[1]]) + abs(y[[2]]))
        })
        return(Reduce(`+`,tmp))
      }
      #minkowski
      multiminkowski <- function(x) {
        tmp <- lapply(x, function(y) {
          abs((y[[1]]-y[[2]])^lambda)
        })
        return(Reduce(`+`,tmp)^(1/lambda))
      }
      #mahalanobis
      multimahalanobis <- function(x){
        tmp <- matrix(unlist(lapply(x,function(y) as.vector(y))),ncol=2)
        tmp <- tmp[!is.na(tmp[,1]),] 
        if( length(tmp)==0 | is.null(dim(tmp)) ) {
          return(NA)
        } else if(rcond(cov(tmp)) <= 0.001) {
          return(NA)
        } else {
          #return the inverse of the covariance matrix of tmp; aka the precision matrix
          inverse<-solve(cov(tmp)) 
          if(debugging){
            print(inverse)
          }
          tmp<-scale(tmp,center=T,scale=F)
          tmp<-as.numeric(t(tmp[1,])%*%inverse%*%tmp[1,])
          return(sqrt(tmp))
        }
      }
      #
      ## Decide what function to use
      #
      if( distance_m=="euclidean") {
        distancef <- get("multieuclidean")
      } else if( distance_m=="manhattan" ) {
        distancef <- get("multimanhattan")
      } else if( distance_m=="canberra" ) {
        distancef <- get("multicanberra")
      } else if( distance_m=="minkowski" ) {
        if( lambda==0 ) {
          stop("The Minkowski Distance for lambda = 0 is Infinity; please choose another value for lambda.")
        } else {
          distancef <- get("multiminkowski") 
        }
      } else if( distance_m=="mahalanobis" ) {
        distancef <- get("multimahalanobis")
        warning("Multimahalanobis distance is not fully supported...")
      }
    } else {
      stop("Distance function not defined for multidimensional Rao's Q; please choose among euclidean, manhattan, canberra, minkowski, mahalanobis!")
    }
    #
    ## Reshape values
    #
    vls<-lapply(input, function(x) {raster::as.matrix(x)})
    #
    ## Rescale and add fake columns and rows for moving w
    #
    hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
    if(rescale) {
      trastersm<-lapply(vls, function(x) {
        t1 <- raster::scale(raster(cbind(ver,rbind(hor,x,hor),ver)))
        t2 <- as.matrix(t1)
        return(t2)
      })
    } else {
      trastersm<-lapply(vls, function(x) {
        cbind(ver,rbind(hor,x,hor),ver)
      })
    }
    #
    ## Loop over all the pixels in the matrices
    #
    if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
      message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, may take some time... \n")
    }
    for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
      for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
        if( length(!which(!trastersm[[1]][c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
          raoqe[rw-w,cl-w] <- NA
        } else {
          tw<-lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
          })
          #
          ## Vectorize the matrices in the list and calculate
          #Among matrix pairwase distances
          lv <- lapply(tw, function(x) {as.vector(t(x))})
          vcomb <- combn(length(lv[[1]]),2)
          vout <- c()
          for(p in 1:ncol(vcomb) ) {
            lpair <- lapply(lv, function(chi) {
              c(chi[vcomb[1,p]],chi[vcomb[2,p]])
            })
            vout[p] <- distancef(lpair)
          }
          raoqe[rw-w,cl-w] <- sum(rep(vout,2) * (1/(window)^4),na.rm=TRUE)
        }
      }
      progress(value=cl, max.value=dim(rasterm)[2]+w, progress.bar = FALSE)
    }
    if(exists("pb")) {
      close(pb) 
      message("\nCalculation of Multidimensional Rao's index complete.\n")
    }
  } else{
    message("Something went wrong when trying to calculate Rao's indiex.")
  }  # end of multimensional RaoQ
  
  #----------------------------------------------------#
  
  #
  ## Shannon
  #
  if( shannon ) {
    message("\nStarting Shannon-Wiener index calculation:\n")
    # Reshape values
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    #
    ## Add "fake" columns and rows for moving window
    #
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
    #
    ## Loop over all the pixels
    #
    for (cl in (1+w):(dim(rasterm)[2]+w)) {
      for(rw in (1+w):(dim(rasterm)[1]+w)) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
          shannond[rw-w,cl-w]<-NA
        } else {
          tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
          if( "NA's"%in%names(tw) ) {
            tw<-tw[-length(tw)]
          }
          tw_values<-as.vector(tw)
          p<-tw_values/sum(tw_values)
          p_log<-log(p)
          shannond[rw-w,cl-w]<-(-(sum(p*p_log)))
        }
      }   
      svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
    } 
    message(("\nCalculation of Shannon's index is also complete!\n"))
  } # End ShannonD
  
  #----------------------------------------------------#
  
  #
  ## Return multiple outputs
  #
  if(debugging){
    message( "#check: return function." )
  }
  
  if( shannon ) {
    if( nc.cores>1 ) {
      outl<-list(do.call(cbind,raop),shannond)
      names(outl)<-c("Rao","Shannon")
      return(outl)
    } else if( nc.cores==1 ){ 
      outl<-list(raoqe,shannond)
      names(outl)<-c("Rao","Shannon")
      return(outl)
    }
  } else if( !shannon & mode=="classic" ) {
    if( isfloat & nc.cores>1 ) {
      return(do.call(cbind,raop)/mfactor)
      if(debugging){
        message("#check: return function - classic.")
      }
    } else if( !isfloat & nc.cores>1 ) {
      outl<-list(do.call(cbind,raop))
      names(outl)<-c("Rao")
      return(outl)
    } else if( isfloat & nc.cores==1 ) {
      outl<-list(raoqe)
      names(outl)<-c("Rao")
      return(outl)    
    } else if( !isfloat & nc.cores==1 ) {
      outl<-list(raoqe)
      names(outl)<-c("Rao")
      return(outl)    
    } else if( !isfloat & nc.cores>1 ) {
      outl<-list(do.call(cbind,raoqe))
      names(outl)<-c("Rao")
      return(outl)
    }
  } else if( !shannon & mode=="multidimension" ) {
    outl<-list(raoqe)
    names(outl)<-c("Multidimension_Rao")
    return(outl)
  }
}
