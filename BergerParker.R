BergerParker <- function(input, window=3, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
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