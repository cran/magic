\name{shift}
\alias{shift}
\alias{ashift}
\title{Shift origin of arrays and vectors}
\description{
Shift origin of arrays and vectors. 
}
  \usage{
shift(x, i)
ashift(a, v)
}
\arguments{
  \item{x}{Vector to be shifted}
  \item{i}{Number of places elements to be shifted}
  \item{a}{Array to be shifted}
  \item{v}{Vector of numbers to be shifted in each dimension}
}
\details{
  Function \code{shift(x,n)} returns \eqn{P^n(x)}{P^n(x)} where
  \eqn{P}{P} is the permutation
  \eqn{(n,1,2,\ldots,n-1)}{(n,1,2,...,n-1)}.

  Function \code{ashift} is the array generalization of this: the
  \eqn{n^{\rm th}}{n-th} dimension is shifted by \code{v[n]}.  In other
  words,
  \code{ashift(a,v)=a[shift(1:(dim(a)[1]),v[1]),...,shift(1:(dim(a)[n]),v[n])]}. 
  It is named by analogy with \code{abind()} and \code{aperm()}.

  This function is here because a shifted semimagic square or hypercube
  is semimagic and a shifted pandiagonal square or hypercube is
  pandiagonal (note that a shifted magic square is not necessarily
  magic, and a shifted perfect hypercube is not necessarily perfect).
}
\author{Robin K. S. Hankin}
\examples{
shift(1:10,3)
m <- matrix(1:100,10,10)
ashift(m,c(1,1))
ashift(m,c(0,1))    #note columns shifted by 1, rows unchanged.
ashift(m,dim(m))    #m unchanged (Mnemonic).
}
\keyword{array}
