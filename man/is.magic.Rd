\name{is.magic}
\alias{is.magic}
\alias{is.panmagic}
\alias{is.semimagic}
\alias{is.associative}
\alias{is.normal}
\title{Various tests for the magicity of a square}
\description{
  Returns TRUE if the square is magic, semimagic, panmagic, associative,
  normal.
}
\usage{
is.magic(m,give.answers=FALSE)
is.panmagic(m,give.answers=FALSE)
is.semimagic(m,give.answers=FALSE)
is.associative(m)
is.normal(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{The square to be tested}
  \item{give.answers}{Boolean, with TRUE meaning return a list of
    rowwise, columnwise and both broken diagonal (NW-SE and NE-SW) sums.}
}
\details{
  A semimagic square is one whose row sums all equal its columnwise
  sums (ie the magic constant).

  A magic square is a semimagic square with the sum of both unbroken
  diagonals equal to the magic constant.

  A panmagic square is a magic square all of whose broken diagonals have
  sum to the magic constant.

  A normal square is one that contains \eqn{n^2}{n^2} consecutive
  integers (typically starting at 0 or 1).  
}
\value{
  Returns TRUE if the square is semimagic, etc.

  If \code{give.answers} is TRUE, return a list of five elements.  The
  first is TRUE if the square is (semimagic, magic, panmagic), elements
  2-5 give the result of a call to \code{allsums()}, viz: rowwise and
  columwise sums; and broken major
  (ie NW-SE) and minor (ie NE-SW) diagonal sums.}
\references{\url{http://mathworld.wolfram.com/MagicSquare.html}}
\author{Robin K. S. Hankin}
\seealso{\code{\link{minmax}},\code{\link{is.perfect}},\code{\link{is.semimagichypercube}}}
\examples{
is.magic(magic(4))
f <- function(n){is.magic(magic(n))}
all(sapply(3:50,f))

is.panmagic(panmagic.4())
is.panmagic(panmagic.8())

proper.magic <- function(m){is.magic(m) & minmax(c(1,diff(sort(m))))}
proper.magic(magic(20))
}
\keyword{array}
