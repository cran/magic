\name{is.semimagichypercube}
\alias{is.semimagichypercube}
\alias{is.magichypercube}
\alias{is.perfect}
\title{magic hypercubes}
\description{
Returns \code{TRUE} if a hypercube is semimagic, magic, perfect
}
\usage{
is.semimagichypercube(a,give.answers=FALSE)
is.magichypercube(a,give.answers=FALSE)
is.perfect(a)
}
\arguments{
  \item{a}{The hypercube (array) to be tested}
  \item{give.answers}{Boolean, with TRUE meaning to also return
    the sums}
}
\details{
  (Although apparently non-standard, here a hypercube is defined to have
  dimension \eqn{d}{d} and order \eqn{n}{n}---and thus has \eqn{n^d}{n^d}
  elements).
  
  A {\em semimagic hypercube} has all ``rook's move'' sums equal to the magic
  constant (that is, each
  \eqn{\sum a[i_1,i_2,\ldots ,i_{r-1},,i_{r+1},\ldots,i_d]}{sum(a[i_1,i_2,\ldots,i_{r-1},,i_{r+1},\ldots,i_d])} with \eqn{1\leq r\leq d}{1 <= r <= d} is
  equal to the magic constant for all values of the \eqn{i}{i}'s).  In
  \code{is.semimagichypercube()}, if \code{give.answers}
  is \code{TRUE}, the sums returned are in the form of an array of
  dimension \code{c(rep(n,d-1),d)}.  The first \code{d-1} dimensions are
  the coordinates of the projection of the summed elements onto the surface
  hypercube.  The last dimension indicates the dimension along which the
  sum was taken over.

  A {\em magic hypercube} is a semimagic hypercube with all
  \eqn{2^{d-1}}{2^(d-1)} major (ie extreme point-to-extreme point)
  diagonals summing correctly.  Note that not all
  subhypercubes are necessarily magic! (for example, consider a
  5-dimensional magic hypercube \code{a}.  The square \code{b} defined
  by \code{a[1,1,1,,]} might not be magic: the diagonals of
  \code{b} are not covered by the definition of a magic hypercube).
  Some subhypercubes of a magic hypercube are not even semimagic: see
  below for an example.

  In \code{is.magichypercube()}, if argument \code{give.answers=TRUE}
  then a list is returned.  The first element of this list is Boolean
  with \code{TRUE} if the array is a magic hypercube.  The second
  element gives the rook's move sums as
  for \code{is.semimagichypercube()} and the third element gives an
  array of \eqn{2^d}{2^d} elements with index 1 meaning the sum is taken
  from elements in \code{a} with the \eqn{d^{\rm th}}{d-th} dimension
  running over \code{1:n} (index 2 means run over \code{n:1}).

  A {\em perfect magic hypercube} is a magic hypercube with all
  nonbroken diagonals summing correctly.  This is a bad-ass
  requirement for high dimensional hypercubes.

  Even in three dimensions, being perfect is pretty bad.  Consider a
  \eqn{5\times5\times 5}{5x5x5} (ie three dimensional), cube.  Say
  \code{a=magiccube.2np1(2)}.  Then the square defined by
  \code{sapply(1:n,function(i){a[,i,6-i]}, simplify=TRUE)}, which is a
  subhypercube of \code{a}, is not even semimagic: the
  rowsums are incorrect (the colsums must sum correctly because \code{a}
  is magic).  Note that the diagonals of this square are two of the
  ``extreme point-to-point'' diagonals of \code{a}.

  A {\em pandiagonal magic hypercube} (or sometimes just a {\em perfect
  hypercube}) has sums equal, and all sums of all diagonals,
  including broken diagonals, are equal.  This is one seriously bad-ass
  requirement.  I reckon that is a total of
  \eqn{\frac{1}{2}\left( 3^d-1\right)\cdot n^{d-1}}{(3^d-1)n^(d-1)/2}
  correct summations.  This is not coded up yet; I can't see how to do
  it in anything like a vectorized manner.

  The terminology in this area is pretty confusing.  
  }
  \references{
    \url{http://members.shaw.ca/hdhcubes/cube_perfect.htm#Cubes}
    "Generalized magic cubes", Dana Richards .  Mathematics
    Magazine, volume 53, number 2, March 1980.
  }
\author{Robin K. S. Hankin}
\seealso{\code{\link{is.magic}}, \code{\link{allsubhypercubes}},\code{\link{hendricks}}}
\examples{
library(abind)
is.semimagichypercube(magiccube.2np1(1))
is.semimagichypercube(magichypercube.4n(1,d=4))

is.perfect(magichypercube.4n(1,d=4))

#Now for a subhypercube of a magic hypercube that is not semimagic:
is.magic(allsubhypercubes(magiccube.2np1(1))[[10]])

data(hendricks)
is.perfect(hendricks)


#note that Hendricks's magic cube also has many broken diagonals summing
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
