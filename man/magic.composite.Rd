\name{magic.composite}
\alias{magic.composite}
\title{Composite magic squares}
\description{
  Gives a magic square that is a product of two magic squares.
}
\usage{
magic.composite(a, b)
}
\arguments{
  \item{a}{First magic square; if \code{a} is an integer, use \code{magic(a)}.}
  \item{b}{as above}
}
\details{
Produces a magic square of order order(a)*order(b).
}
\references{
  William H. Benson and Oswald Jacoby. New recreations with magic squares, Dover 1976.}
\author{Robin K. S. Hankin}
\examples{
magic.composite(3,4)
magic.composite(4,3)
}
\keyword{array}
