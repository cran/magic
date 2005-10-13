\name{magic-package}
\alias{magic-package}
\docType{package}
\title{
Magic squares and hypercubes; arbitrary dimensioned array manipulation
}
\description{

A collection of efficient, vectorized algorithms for the creation and
investigation of magic squares and hypercubes, including a variety of
functions for the manipulation and analysis of arbitrarily dimensioned
arrays.

The package includes methods for creating normal magic squares of any
order greater than 2.  The ultimate intention is for the package to be a
computerized embodiment all magic square knowledge, including direct
numerical verification of properties of magic squares (such as recent
results on the determinant of odd-ordered semimagic squares).

}
\details{
\tabular{ll}{
Package: \tab magic\cr
Type: \tab Package\cr
Version: \tab 1.3-21\cr
Date: \tab 2005-11-30\cr
License: \tab GPL\cr
}
}
\author{
Robin K. S. Hankin

Maintainer: r.hankin@noc.soton.ac.uk
}
\keyword{ package }
\examples{
magic(5)

a <-  magiccube.2np1(1)
adiag(1,a)
apad(a,2,1)
allsubhypercubes(a)
arev(a)
apltake(a,c(2,2))
arot(a)
arow(a,1)
}
