\name{magic}
\alias{magic}
\title{Creates magic squares}
\description{
  Creates magic squares of any order \eqn{>2}{>2}.  Uses the appropriate
method depending on n modulo 4.
}
\usage{
magic(n)
}
\arguments{
  \item{n}{Order of magic square}
}
\details{
Calls either \code{magic.2np1()}, \code{magic.4n()},
or \code{magic.4np2()} depending on the value of \code{n}.  Returns a
magic square in standard format.
}
\references{William H. Benson and Oswald Jacoby.  New recreations with
  magic squares, Dover 1976.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{magic.2np1}}, \code{\link{magic.prime}},
  \code{\link{magic.4np2}}, \code{\link{magic.4n}},
  \code{\link{as.standard}}, \code{\link{force.integer}}}
\examples{
magic(6)
all(sapply(3:10,function(n){is.magic(magic(n))}))
}
\keyword{array}
