\name{ChIPsim-package}
\Rdversion{1.1}
\alias{ChIPsim-package}
\alias{ChIPsim}
\docType{package}
\title{
Simulation of ChIP-seq experiments
}
\description{
This package provides a framework for the simulation of ChIP-seq experiments. An implementation of
a simulation for nucleosome positioning experiments is provided as part of the package. Simulations
for other experiments can be implemented using the provided framework.
}
\details{
\tabular{ll}{
Package: \tab ChIPsim\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2009-09-20\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
Depends: \tab Biostrings\cr
Imports: \tab IRanges, ShortRead\cr
Suggests: \tab actuar, zoo\cr
}
Function \code{\link{simChIP}} is the main driver of the simulation. To simulate different
types of experiments the functions passed to the \code{functions} argument of \code{\link{simChIP}}
have to be replaced. See the vignettes for more detail.
}
\author{
Peter Humburg

Maintainer: Peter Humburg <Peter.Humburg@csiro.au>
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
\code{\link[ShortRead:ShortRead-package]{ShortRead}} and its dependencies 
are used to handle short read and genomic sequences.
}
\examples{
## See the accompanying vignette 'Introduction to ChIPsim' for a detailed 
## example of how to use this package for nucleosome positioning simulations. 
## A guide for the writing of extensions that cover other types of
## experiments is provided in 'Extending ChIPsim'.
}
