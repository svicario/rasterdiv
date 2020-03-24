\name{CRE}
\alias{CRE}
\title{Cumulative Residual Entropy  (CRE)}
\description{
Computes Cumulative Residual Entropy  (CRE) on different classes of numeric matrices using a moving window algorithm.
}
\usage{
CRE(x, window=9, mode="classic",  
shannon=FALSE, rescale=FALSE, na.tolerance=0.0, simplify=3, 
n.processes=1, cluster.type="SOCK", debugging=FALSE)
}
\arguments{
  \item{x}{input data may be a matrix, a Spatial Grid Data Frame, a RasterLayer or a list of these objects. In the latter case, if \emph{mode="classic"} only the first element of the list will be considered.}
  \item{window}{the side of the square moving window, it must be a odd numeric value greater than 1 to ensure that the target pixel is in the center of the moving window. Default value is 3. If \emph{proxy::dist} is a matrix then the function will assume that this is the distance matrix, and therefore no distance will be derived.}
  \item{mode}{currently, there are two modes to calculate Cumulative Residual Entropy (CRE). If mode is "classic", then the monodimension CRE will be calculated on one single matrix. If mode is "multidimension" (experimental!) a list of matrices must be provided as input. In this latter case, the multidimensional cumulative residual probability will be calculated . Default value is "classic"}
    \item{shannon}{a boolean variable set to FALSE by default. If TRUE, a matrix with Shannon index will be also calculated.  Default value is FALSE.}
  \item{rescale}{a boolean variable set to FALSE by default. If TRUE, \emph{x} will be scaled and centered to standardise different matrices if mode is "multidimension". Default value is FALSE.}
  \item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate CRE in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise CRE will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
  \item{simplify}{Number of decimal places to be retained to calculate CRE. Only if \emph{x} is floats.}
  \item{n.processes}{the number of processes (cores) which will be spawned. Default value is 1.}
  \item{cluster.type}{the type of cluster which will be created. The options are \emph{"MPI"} (which calls "makeMPIcluster"), \emph{"FORK"} and \emph{"SOCK"} (which call "makeCluster"). Default type is \emph{"SOCK"}.}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
  }
\details{
Unidimension Cumulative Residual Entropy (\eqn{CRE}) is calculated on a numerical vector  as \eqn{CRE = CRE=-\sum_{i=1}^{N}{P(X\geq x_i)\log{P(X\geq x_i)dx}} [1] where  \eqn{ dx=(x_i -x_{i-1})} and \emph{P} is  vector the probability that values of the vector of observation are larger or equal to each value of the vector. 
	 In the "multidimension" CRE, \emph{P} becomes an array with as dimension as the one of the observation. In each cell of \emph{P} the probability, estimates with the frequency of the number of observation that as the same time satisify the larger or equal requirment for the a different combination of values along the dimension. dx becomes the products of the difference along each dimension. \cr
The theoretical minimum of CRE is 0, when all values are identical in a window.  The values of CRE increase with the range of the observations and the more the observation are equally spread across values.}
\value{
If shannon=TRUE, a list of matrices of length two, otherwise a matrix of dimension \emph{dim(x)}.
}
\references{
[1] Rao M, Chen Y, Vemuri BC, Wang F (2004) Cumulative Residual Entropy: A New Measure of Information. IEEE Trans Inf Theory 50:1220–1228.
}
\author{
saverio vicario \email{saverio.vicario@cnr.it}\cr
}

\examples{
#Minimum example; compute Rao's index 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
out <- CRE(x=a,window=3,na.tolerance=0,shannon=FALSE,mode="classic")

#Compute both Rao and Shannon index 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
out <- CRE(x=a,window=3,na.tolerance=0,shannon=TRUE,mode="classic")

#Compute multidimension Rao's index rescaling the input matrices
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
b <- matrix(c(0.5,0.5,0.1,0.1,0.3,0.3,0.3,0.3,0.3),ncol=3,nrow=3)
out <- CRE(x=list(a,b),window=3,na.tolerance=0,
  mode="multidimension",rescale=TRUE,debugging=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{Cumulative Residual Entropy}% __ONLY ONE__ keyword per line
\keyword{CRE}% __ONLY ONE__ keyword per line
