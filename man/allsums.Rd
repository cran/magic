\name{allsums}
\alias{allsums}
\title{Row, column, and two diagonal sums of arrays}
\description{
  Returns all rowsums, all columnsums, and all (broken) diagonal sums of
  a putative magic square.
}
\usage{
allsums(m)
}
\arguments{
  \item{m}{The square to be tested}
}
\value{
  Returns a list of four elements.
  \item{rowsums}{All n row sums}
  \item{colsums}{All n column sums}
  \item{majors}{All n broken major diagonals (northwest-southeast).
    First element is the long (unbroken) major diagonal, tested by
    \code{is.magic()}}
  \item{minors}{All n broken minor diagonals (northeast-southwest).
    First element is the long (unbroken) minor diagonal.}
}
\author{Robin K. S. Hankin}
\note{Used by \code{is.magic()} et seq}
\seealso{\code{\link{is.magic}},\
\code{\link{is.semimagic}},\
\code{\link{is.panmagic}}
  }
\examples{
allsums(magic.prime(7))
}
\keyword{array}
