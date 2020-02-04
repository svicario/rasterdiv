\name{spectralrao}
\alias{spectralrao}
\title{
Rao computation function
}
\description{
Code to calculate Rao's quadratic entropy on a numeric matrix, RasterLayer object (or lists) using a moving window algorithm. The function also calculates Shannon-Wiener index.
}
\usage{
spectralrao(input, distance_m="euclidean", p=NULL, window=9, mode="classic", lambda=0, 
shannon=FALSE,rescale=FALSE, na.tolerance=0.0, simplify=3, 
nc.cores=1, cluster.type="MPI", debugging=FALSE, ...)
}
\arguments{
  \item{input}{the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.}
  \item{distance_m}{}
  \item{p}{}
  \item{window}{the size of the square window size. Default value is 9.}
  \item{mode}{}
  \item{lambda}{}
  \item{shannon}{}
  \item{rescale}{}
  \item{na.tolerance}{}
  \item{simplify}{the 10 power which will be used to covert float into natural number. Default value is 3.}
  \item{nc.cores}{the nuber of cores which will be used. Default value is 1.}
  \item{cluster.type}{the type of cluster which will be used. Default type is "MPI".}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, let the user check all the steps}
}
\details{
Rao's Q Min = 0, if all pixel classes have distance 0. If the chosen distance ranges between 0 and 1, Rao's Max = 1-1/S (Simpson Diversity, where S is the number of pixel classes).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matteo Marcantonio \email{}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
