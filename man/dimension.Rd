\name{dimension}
\alias{dimension}
\title{Generalized row and col}
\description{
Given an array, returns an array of the same size whose elements
are sequentially numbered along the \eqn{i^{\rm th}}{i-th} dimension.

Used by \code{magiccube.2np1}.
}
\usage{
dimension(a, i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{array to be converted}
  \item{i}{Number of the dimension}
}
\value{
An integer matrix with the same dimensions as \code{a}, with element 
\eqn{\left(n_1,n_2,\ldots n_d\right)^{\rm th}}{(n1,n2, ..., n_d)}
being \eqn{n_i}{n_i}.
}
\author{Robin K. S. Hankin}
\examples{
a <- array(0,c(3,3,2,2))
dimension(a,2)
a[(dimension(a,1)+dimension(a,2)+dimension(a,3)+dimension(a,4))\%\%2==1] <- 1
a
}
\keyword{array}
