\name{Hill}
\alias{Hill}
\title{Hill's index of diversity - Hill numbers (D)}
\description{
Computes Hill's index of diversity (Hill numbers) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
Hill(x, window = 3, alpha = 1, n.processes = 1, cluster.type = "SOCK", debugging = FALSE)
  }
  \arguments{
  \item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, only the first element of the list will be considered.}
  \item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3.}
  \item{alpha}{Order of the Hill number to compute the index. If \emph{"alpha"} is a vector with length greater than 1, then the index will be calculated over \emph{x} for each value in the sequence.}
  \item{n.processes}{the number of processes (cores) which will be spawned. Default value is 1.}
  \item{cluster.type}{the type of cluster which will be created. The options are "MPI" (calls "makeMPIcluster"), "FORK" and "SOCK" (call "makeCluster"). Default type is "SOCK".}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
  }
  \details{
Hill numbers (\eqn{{}^qD}) are calculated on a numerical matrices as \eqn{{}^qD = (\sum_{i=1}^{R} {p^q}_i)^{1/(1-q)}},where \emph{q} is the order of the Hill number, \emph{R} is the total number of categories (i.e., unique numerical values in a numerical matrix), \emph{p} is the relative abundance of each category. When q=1, Shannon.R is called to calculate \eqn{exp(H^1)} instead of the indefinite \eqn{{}^1D}. if \eqn{q > 2*10^9}, BerkgerParker.R is called to calculate \eqn{1/{{}^\infty D}}. Hill numbers of low order weight more rare categories, whereas Hill numbers of higher order weight more dominant categories.
}
\value{
A list of matrices with length equal to the length of "alpha".
}
\references{
Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology 54, 427-431.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\note{
Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi

Microsoft Windows users may need some additional work to use "MPI", see:\cr
https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as , ~~~
\code{\link{BergerParker}}
\code{\link{Shannon}}
}
\examples{
#Minimum example; compute Hill's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
hill <- Hill(x=a,window=3,alpha=1:5)

#Compute Hill's index with alpha=2 on global Copernicus NDVI using parallel computing
%not running%
# require("proxy")
# require("raster")
# require("svMisc")

# ndvi17_r <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE)
# hill <- Hill(x=ndvi17_r,window=9,n.processes=4,cluster.type="SOCK",alpha=0:5)
# hillra <- lapply(hill, function(x) raster(x,template=ndvi17))
# cols <- colorRampPalette(brewer.pal(9,"RdBu"))

#Plot results using lattice

# levelplot(hillra,main="Hill's number (9px moving window) - Copernicus NDVI 5km",
#  col.regions = cols,as.table = T,layout = c(0,5,1),
#  names.attr = paste("alpha",0:5,sep=" "))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods.}% use one of  RShowDoc("KEYWORDS")
\keyword{Hill's index.}% __ONLY ONE__ keyword per line
\keyword{Hill's numbers.}% __ONLY ONE__ keyword per line
\keyword{Index of diversity.}% __ONLY ONE__ keyword per line
