\name{magic.constant}
\alias{magic.constant}
\title{Magic constant of a magic square}
\description{
  Returns the magic consant: that is, the common sum for all rows,
  columns and (broken) diagonals of a magic square
}
\usage{
magic.constant(n)
}
\arguments{
  \item{n}{Order of the square}
}
\details{
  Returns \eqn{n(n^2+1)/2}{n(n^2+1)/2}.  Note that this is correct for
  magic squares whose numbers run from \eqn{1\ldots n^2}{1...n^2}.  If the
  square is numbered \eqn{0\ldots n^2-1}{1 ... n^2-1} the correct magic
  constant would be \eqn{n(n^2-1)/2}{n(n^2-1)/2}, or
  \code{magic.constant(n)-n}. 
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{magic}}}
\examples{
magic.constant(4)
}
\keyword{array}
