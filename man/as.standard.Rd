\name{as.standard}
\alias{as.standard}
\alias{is.standard}

\title{Standard form for magic squares}
\description{
  Transforms a magic square or magic hypercube into Frenicle's standard form
}
\usage{
as.standard(a)
is.standard(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Magic square or hypercube (array) to be tested or transformed}
}

\details{
\code{as.standard()} transforms a magic square into Frenicle's standard
form.  The four numbers at each of the four corners are determined.
First, the square is rotated so the smallest of the four is at the upper
left.  Then, element \code{[1,2]} is compared with element\code{[2,1]}
and, if it is larger, the transpose is taken.

Thus all eight rotated and transposed versions of a magic square have
the same standard form.

The square returned by \code{magic()} is in standard form.


For hypercubes, the algorithm is generalized.  Firstly, the hypercube is
reflected so that \code{a[1,1, ... ,1,1]} is the smallest of the
\eqn{2^d}{2^d} corner elements (eg \code{a[1,n,1, ... ,1,1]}).

Next, \code{aperm()} is called so that \code{a[1,1,  ... ,1,2] < a[1,1, ... ,2,1] < ... <  a[2,1, ... ,1,1]}
  
\code{is.standard()} returns \code{TRUE} if the magic square or
hypercube is in standard form.  \code{is.standard()} and
\code{as.standard()} check for neither magicness nor normality (use
\code{\link{is.magic}} and \code{\link{is.normal}} for this).
}

\author{Robin K. S. Hankin}
\seealso{\code{\link{magic}},\ \code{\link{eq}}}
\examples{
is.standard(magic.2np1(4))
as.standard(magic.4n(3))

as.standard(magichypercube.4n(1,5))

}
\keyword{array}
