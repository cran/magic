\name{allsums}
\alias{allsums}
\title{Row, column, and two diagonal sums of arrays}
\description{
  Returns all rowsums, all columnsums, and all (broken) diagonal sums of
  a putative magic square.
}
\usage{
allsums(m,FUN=sum)
}
\arguments{
  \item{m}{The square to be tested}
  \item{FUN}{Function, defaulting to \code{sum}, to be applied to the
  square rowwise, columnwise, and diagonalwise}
}
\value{
  Returns a list of four elements.  In the following, \dQuote{sums}
  means \dQuote{the result of applying FUN()}.
  \item{rowsums}{All n row sums}
  \item{colsums}{All n column sums}
  \item{majors}{All n broken major diagonals (northwest-southeast).
    First element is the long (unbroken) major diagonal, tested by
    \code{is.magic()}}
  \item{minors}{All n broken minor diagonals (northeast-southwest).
    First element is the long (unbroken) minor diagonal.}
}
\author{Robin K. S. Hankin}
\note{
If \code{FUN()} returns a vector, then the columns of the four elements
are the result of applying \code{FUN()}.  See third example below.

Used by \code{is.magic()} et seq.
}
\seealso{\code{\link{is.magic}},\
\code{\link{is.semimagic}},\
\code{\link{is.panmagic}}
  }
\examples{
allsums(magic(7))
allsums(magic(7),FUN=max)
allsums(magic(7),FUN=function(x){x[1:2]})
  # shows how the minor diagonals are ordered: first [1,n] to [n,1] then
  # [1,n+1] to [n,2] etc.

allsums(magic.prime(7),sort)
  # beware! compare apply(magic(7),1,sort) with apply(magic(7),2,sort)



}
\keyword{array}
