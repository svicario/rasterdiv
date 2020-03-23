\name{Rao}
\alias{Rao}
\title{Rao's index of quadrtic entropy (Q)}
\description{
Computes Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
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
  \item{na.tolerance}{a numeric value \eqn{(0.0-1.0)} which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default values is 0.0 (i.e., no tolerance for NA's).}
  \item{simplify}{Number of decimal places to be retained to calculate distances in Rao's index. Only if \emph{x} is floats.}
  \item{n.processes}{the number of processes (cores) which will be spawned. Default value is 1.}
  \item{cluster.type}{the type of cluster which will be created. The options are \emph{"MPI"} (which calls "makeMPIcluster"), \emph{"FORK"} and \emph{"SOCK"} (which call "makeCluster"). Default type is \emph{"SOCK"}.}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.}
  }
\details{
Classical Rao's Index (\eqn{Q}) is calculated on a numerical matrix as \eqn{Q = \sum_{i=1}^{R}{\sum_{j=1}^{R}{d_{i,j}\times{p_i}\times{p_j}}}} [1]. Where R is the number of categories, whereas \emph{i} and \emph{j} are pair of numerical categories in the same moving window. In the "multidimension" Rao's index, distances among categories are calculated considering more than one layer, then the pairwise distance between each pair of numerical categories is multiplied to the square of the size of the moving window (which is somewhat the same as to calculate the variance of the multidimensional distance [2].). \cr
The theoretical minimum of Rao's Q is 0, when all categories in a window have distance 0. If the distance chosen to calculate Rao's Index ranges between 0 and 1, the maximum value of Rao's Index equals the Simpson Index of Diversity \eqn{1-1/S_i} where \emph{S} is the number of categories in window \emph{i}.}
\value{
If shannon=TRUE, a list of matrices of length two, otherwise a matrix of deimension \emph{dim(x)}.
}
\references{
[1] Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 2443.
[2] Rocchini, D., M. Marcantonio, and C. Ricotta (2017). Measuring Rao’s Q diversity index from remote sensing: An open source solution. Ecological Indicators. 72: 234–238.
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
