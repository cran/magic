\name{is.magic}
\alias{is.magic}
\alias{is.panmagic}
\alias{is.semimagic}
\alias{is.associative}
\alias{is.normal}
\alias{is.mostperfect}
\alias{is.2x2.correct}
\alias{is.diagonally.correct}
  
\title{Various tests for the magicness of a square}
\description{
  Returns TRUE if the square is magic, semimagic, panmagic, associative,
  normal.  If argument \code{give.answers} is TRUE, also returns
  additional information about the sums.
}
\usage{
is.magic(m,give.answers=FALSE)
is.panmagic(m,give.answers=FALSE)
is.semimagic(m,give.answers=FALSE)
is.associative(m)
is.normal(m)
is.mostperfect(m,give.answers=FALSE)
is.2x2.correct(m,give.answers=FALSE)
is.diagonally.correct(m,give.answers=FALSE)
}
\arguments{
  \item{m}{The square to be tested}
  \item{give.answers}{Boolean, with TRUE meaning return additional
    information about the sums (see details).}
}
\details{
  A {\em semimagic square} is one all of whose row sums equal all its
  columnwise sums (ie the magic constant).

  A {\em magic square} is a semimagic square with the sum of both unbroken
  diagonals equal to the magic constant.

  A {\em panmagic square} is a magic square all of whose broken diagonals
  sum to the magic constant.  Ollerenshaw calls this a ``pandiagonal''
  square.

  A {\em most perfect} square has all 2-by-2 arrays anywhere within the square
  summing to \eqn{2S} where \eqn{S=n^2+1}{S=n^2+1}; and all pairs of
  integers \eqn{n/2}{n/2} distant along the same major (NW-SE) diagonal sum to
  \eqn{S}{S} (note that the \eqn{S}{S} used here differs
  from Ollerenshaw's because her squares are numbered starting at zero).
  The first condition is tested by  \code{\link{is.2x2.correct}} and the
  second by \code{\link{is.diagonally.correct}}.
  
  All most perfect squares are panmagic.

  A {\em normal square} is one that contains \eqn{n^2}{n^2} consecutive
  integers (typically starting at 0 or 1).  
}
\value{
  Returns TRUE if the square is semimagic, etc.
 
  If \code{give.answers} is taken as an argument and is TRUE, return a
  list of at least five elements.  The first element of the list is the
  answer: it is TRUE if the square is (semimagic, magic, panmagic).
  Elements 2-5 give the result of a call to \code{allsums()}, viz: rowwise
  and columwise sums; and broken major (ie NW-SE) and minor (ie NE-SW)
  diagonal sums.

  Function \code{is.diagonally.correct()} also returns the sums of
  elements distant \eqn{n/2}{n/2} along a major diagonal
  (\code{diag.sums}); and function \code{is.2x2.correct()} returns the
  sum of each \eqn{2\times 2}{2x2} submatrix (\code{tbt.sums}).
  Function \code{is.mostperfect()} returns both of these.


}
\references{\url{http://mathworld.wolfram.com/MagicSquare.html}}
\author{Robin K. S. Hankin}
\seealso{\code{\link{minmax}},\code{\link{is.perfect}},\code{\link{is.semimagichypercube}}}
\examples{
is.magic(magic(4))
f <- function(n){is.magic(magic(n))}
all(sapply(3:50,f))

is.panmagic(panmagic.4())
is.panmagic(panmagic.8())

data(Ollerenshaw)
is.mostperfect(Ollerenshaw)


proper.magic <- function(m){is.magic(m) & minmax(c(1,diff(sort(m))))}
proper.magic(magic(20))
}
\keyword{array}
