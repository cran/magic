\name{allsubhypercubes}
\alias{allsubhypercubes}
\title{Subhypercubes of magic hypercubes}
\description{
Extracts all subhypercubes from an n-dimensional hypercube.
}
\usage{
allsubhypercubes(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{The magic hypercube whose subhypercubes are computed}
}
\value{
  Returns a list, each element of which is a subhypercube.  Note that
  major diagonals are also returned (as n-by-1 arrays).
}
\author{Robin K. S. Hankin}
\note{This function is a dog's dinner.  It's complicated, convoluted,
  and needs two absurd uses of the \code{eval(parse(text=...))}
  construction.  Basically it sucks big time.

  BUT\ldots I cannot for the life of me see a better way that gives the
  same results, without loops, on hypercubes of arbitrary dimension.

  The function cannot cope with \code{d} as high as 5, for \code{n=4}.
  Heigh ho.
  }
\seealso{\code{\link{is.perfect}}}
\examples{
  library(abind)
  a <- magichypercube.4n(1,d=4)
  allsubhypercubes(a)
}
\keyword{array}
