Hill <- function(input, window=3, mode ="single", BergerParker=FALSE, alpha=1, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...){
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