\name{recurse}
\alias{recurse}
\title{Recursively apply a permutation}
\description{
Recursively apply a permutation to a vector an arbitrary number of
times.  Negative times mean apply the inverse permutation.
}
\usage{
recurse(start = 1:n, perm, i)
}
\arguments{
  \item{start}{Start vector to be permuted}
  \item{perm}{Permutation (integers 1 to \code{length(start)} in some order)}
  \item{i}{Number of times to apply the permutation.  \code{i=0} gives
  \code{start} by definition}
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hudson}}}
\examples{
recurse(start=rnorm(10),perm=sample(10),i=3)
}
\keyword{array}
