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