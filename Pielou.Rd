\name{Pielou}
\alias{Pielou}
\title{
Pielou computation function
}
\description{
Compute the Pielou index, after data check
}
\usage{
Pielou(input, window=3, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE,   ...)
}
\arguments{
 \item{input}{the function input could be a matrix, a Spatial Grid data frame, a raster layer or list of matrix. The index will be compted over it.}
 \item{window}{the size of the square window size. Default value is 3.}
 \item{simplify}{the 10 power which will be used to covert float into natural number. Default value is 3.}
 \item{nc.cores}{the nuber of cores which will be used. Default value is 1.}
 \item{cluster.type}{the type of cluster which will be used. Default type is "MPI".}
 \item{debugging}{a boolean variable set to FALSE by default. If TRUE, let the user check all the steps}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix or a list of matrixes with the Pielou index computed through moving window of the given size
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
