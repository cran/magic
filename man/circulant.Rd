\name{circulant}
\alias{circulant}
\alias{is.circulant}
\title{Circulant matrices of any order}
\description{
Circulant matrices of any order
}
\usage{
circulant(n, vec=1:n)
is.circulant(m)
}
\arguments{
  \item{n}{Order of circulant matrix}
  \item{vec}{Vector of elements of the first row, defaulting to
    \code{1:n}.}
  \item{m}{Matrix to be tested for circulantism}
}
\details{
  A matrix \eqn{a}{a} is circulant if all (major) diagonals are uniform,
  ie if  \eqn{a_{ij}=a_{kl}}{a[i,j]==a[k,j]} when \eqn{i-j=k-l}{i-j=k-l}
  (modulo \eqn{n}{n}).  The standard values to use give \code{1:n}
  for the top row.
}
\references{Arthur T. Benjamin and K. Yasuda.  \emph{Magic
    \dQuote{Squares} Indeed!}, American Mathematical Monthly, vol
  106(2), pp152-156, Feb 1999}
\author{Robin K. S. Hankin}
\examples{
circulant(5)
circulant(5,vec=2^(0:4))
is.circulant(circulant(5))
}
\keyword{array}
