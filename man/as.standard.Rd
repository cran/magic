\name{as.standard}
\alias{as.standard}
\alias{is.standard}

\title{Standard form for magic squares}
\description{
  Transforms a magic square into Frenicle's standard form.
}
\usage{
as.standard(m)
is.standard(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Magic square to be tested or transformed}
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

\code{is.standard()} returns \code{TRUE} if the magic square is in
standard form.  This function checks for neither magicness nor normality
(use \code{\link{is.magic}} and \code{\link{is.normal}} for this).
}

\author{Robin K. S. Hankin}
\seealso{\code{\link{magic}},\ \code{\link{eq}}}
\examples{
is.standard(magic.2np1(4))
as.standard(magic.4n(3))
}
\keyword{array}
