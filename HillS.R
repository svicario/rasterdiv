HillS <- function(rasterm, w, alpha, debugging){
  message("\nStarting Hill index calculation with parameter value = ",alpha," \n\n\n\n")
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