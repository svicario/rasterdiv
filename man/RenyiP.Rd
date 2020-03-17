\name{RenyiP}
\alias{RenyiP}
\title{
Parallelised Renyi's diversity index.
}
\description{
Parallelised Renyi's diversity index.
}
\usage{
RenyiP(rasterm, w, alpha, base, debugging)
}
\arguments{
\item{rasterm}{input data.}
\item{w}{half of the side of the square moving window.}
\item{alpha}{alpha value for order of diversity in Renyi's Index.}
\item{base}{base of the logarigthm to calculate the index.}
\item{debugging}{a boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.}
}
\value{
Matrix or a list of matrixes with the Renyi index computed through moving window of the given size.
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Martina Iannacito \email{martina.iannacito@inria.fr} \cr
Duccio Rocchini \email{duccio.rocchini@unibo.it} \cr
}
\seealso{
\code{\link{Renyi}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
