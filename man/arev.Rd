\name{arev}
\alias{arev}
\title{Reverses some dimensions; a generalization of rev}
\description{
A multidimensional generalization of \code{rev()}: given an array
\code{a}, and a Boolean vector \code{swap}, return an array of the
same shape as \code{a} but with dimensions corresponding to \code{TRUE}
elements of \code{swap} reversed.
}
\usage{
arev(a, swap = TRUE)
}
\arguments{
  \item{a}{Array to be reversed}
  \item{swap}{Vector of Boolean variables.  If \code{swap[i]} is
    \code{TRUE}, then dimension \code{i} of array \code{a} is reversed.
    If \code{swap} is of length one, recycle to \code{length(dim(a))}.}
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{ashift}}}
\examples{
a <- matrix(1:42,6,7)
arev(a)  #Note swap defaults to TRUE

b <- magichypercube.4n(1,d=4)
arev(b,c(TRUE,FALSE,TRUE,FALSE))
}
\keyword{array}
