\name{magic.4n}
\alias{magic.4n}
\title{Magic squares of order 4n}% \eqn{4n}{4n}}
\description{
Produces an associative magic square of order \eqn{4n}{4n} using the delta-x method.
}
\usage{
magic.4n(m, method="a")
}
\arguments{
  \item{m}{Order \eqn{n}{n} of the square is given by \eqn{n=4m}{n=4m}}
  \item{method}{string specifying which method to use.  Method ``a''
    appears to be uniformly better, and clearer.}
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{magic}}}
\examples{
magic.4n(4)
is.magic(magic.4n(5))
}
\keyword{array}

