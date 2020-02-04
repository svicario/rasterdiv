\name{Renyi}
\alias{Renyi}
\title{
Renyi computation function
}
\description{
Compute the Renyi index, after data check
}
\usage{
Renyi(input, window=3, mode ="single", BergerParker=FALSE, alpha=1, base=exp(1), simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...)
}
\arguments{
  \item{input}{the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.}
  \item{window}{the size of the square window size. Default value is 3.}
  \item{mode}{a variable which can be \emph{"single"} to compute Renyi index for just one alpha value, the dealt value; \emph{"iterative"} to compute Renyi index for all the integers value of alpha in a given interval, or \emph{"sequential"} to compute Renyi index for all the alpha values in a given vector.}
  \item{BergerParker}{a boolean variable. The default value is FALSE. If TRUE, in the "single" mode the function will return just the Berger-Parker index. If TRUE, in all the other modes, it will return also the Berger-Parker index.}
  \item{simplify}{the 10 power which will be used to covert float into natural number. Default value is 3.}
  \item{alpha}{a numerical variable. Its default value is 1. In \emph{"single"} mode, alpha has to be a numerical value greater than 0; in \emph{"iterative"} mode, alpha has to be a length 2 vector and in \emph{"sequential"} mode, alpha has to be a vector of leangth at least 2.}
  \item{base}{a numerical value, which lets the user choose the base of the logarithm in R?nyi index formula. Its defalut value is exp(1)} 
  \item{nc.cores}{the nuber of cores which will be used. Default value is 1.}
  \item{cluster.type}{the type of cluster which will be used. Default type is "MPI".}
  \item{debugging}{a boolean variable set to FALSE by default. If TRUE, let the user check all the steps}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix or a list of matrixes with the Renyi index computed through moving window of the given size
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
