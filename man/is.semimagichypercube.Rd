\name{is.semimagichypercube}
\alias{is.semimagichypercube}
\alias{is.magichypercube}
\alias{is.perfect}
\title{magic hypercubes}
\description{
Returns \code{TRUE} if a hypercube is semimagic.
}
\usage{
is.semimagichypercube(a)
is.magichypercube(a)
is.perfect(a)
}

\arguments{
  \item{a}{The hypercube (array) to be tested.}
}
\details{

  (Although apparently non-standard, here a hypercube is defined to have
  dimension \eqn{d}{d} and order \eqn{n}{n}---and thus has \eqn{n^d}{n^d}
  elements).
  
  A semimagic hypercube has all ``rook's move'' sums equal to the magic
  constant.

  A magic hypercube is a semimagic hypercube with all
  \eqn{2^{d-1}}{2^(d-1)} major (ie extreme point-to-extreme point)
  diagonals summing correctly.  Note that not all
  hyperplane subsets are necessarily magic! (for example, consider a
  5-dimensional magic hypercube \code{a}.  The square \code{b} defined
  by \code{a[1,1,1,,]} might not be semimagic: the rows and columns of
  \code{b} are not covered by the definition of a magic hypercube.

  A perfect magic hypercube is a magic hypercube with all nonbroken
  diagonals summing correctly.  This is a bad-ass requirement for high
  dimensional hypercubes.

  Even in three dimensions, being perfect is pretty bad.  Consider a
  \eqn{5\times5\times 5}{5x5x5} (ie three dimensional), cube.  Say
  \code{a=magiccube.2np1(2)}.  Then the square defined by
  \code{sapply(1:n,function(i){a[,i,6-i]}, simplify=TRUE)}, which is a
  subhypercube of \code{a}, is not even semimagic: the
  rowsums are incorrect (the colsums must sum correctly because \code{a}
  is magic).  Note that the diagonals of this square are two of the
  ``extreme point-to-point'' diagonals of \code{a}.

  A pandiagonal magic hypercube (or sometimes just a perfect hypercube)
  has sums equal, and all sums of all diagonals, including broken
  diagonals, are equal.  This is one seriously bad-ass requirement.  I
  reckon that is a total of
  \eqn{\frac{1}{2}\left( 3^d-1\right)\cdot n^{d-1}}{(3^d-1)n^(d-1)/2}
  correct summations.  I can't see how to code this up in anything like
  a vectorized manner.

  The terminology in this area is pretty confusing.  
  }
\references{\code{http://members.shaw.ca/hdhcubes/cube_perfect.htm#Cubes}}
\author{Robin K. S. Hankin}
\seealso{\code{\link{is.magic}}, \code{\link{allsubhypercubes}},\code{\link{hendricks}}}
\examples{
library(abind)
is.semimagichypercube(magiccube.2np1(1))
is.semimagichypercube(magichypercube.4n(1,d=4))

is.perfect(magichypercube.4n(1,d=4))

data(hendricks)
is.perfect(hendricks)


#note that Hendricks' magic cube also has many broken diagonals summing
#correctly:

a <- allsubhypercubes(hendricks)
ld <- function(a){length(dim(a))}

jj <- unlist(lapply(a,ld))
f <- function(i){is.perfect(a[[which(jj==2)[i]]])}
all(sapply(1:sum(jj==2),f))

#but this is NOT enough to ensure that it is pandiagonal (but I
#think hendricks is pandiagonal).

}
\keyword{array}
