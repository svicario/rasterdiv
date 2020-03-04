BergerParker <- function(input, window=3, n.process=1, cluster.type="MPI", debugging=FALSE,   ...){
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
    rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
  }
  #Print user messages
  else if( is(input,"matrix") | is(input,"list") ) {
    rasterm<-as.matrix(rasterm)
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
  if (n.process == 1){
    outS <- BergerParkerS(rasterm, w, debugging)
    message(("\nCalculation of Berger-Parker's index is complete!\n"))
    return (outS)
  }
  else {
    message("##################### Starting parallel calculation #######################")
    #
    if(debugging){cat("#check: Berger-Parker parallel function.")}
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(n.process,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(n.process,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
    if(isfloat) {
      parallel::clusterExport(cl=cls, varlist=c("mfactor"))
    }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()
    outP <- do.call(cbind,BergerParkerP(rasterm, w,  debugging))
    return(outP)
  }
}