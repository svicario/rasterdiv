\name{ShannonS}
\alias{ShannonS}
\title{
Shannon's diversity index.
}
\description{
Shannon's diversity index.
}
\usage{
ShannonS(rasterm, w, debugging)
}
\arguments{
\item{rasterm}{input data.}
\item{w}{half of the side of the square moving window.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.}
}
\value{
Matrix or a list of matrixes with the Shannon index computed through moving window of the given size.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\seealso{
\code{\link{Shannon}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
