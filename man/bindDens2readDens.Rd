\name{bindDens2readDens}
\Rdversion{1.1}
\alias{bindDens2readDens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a feature density into a read density
}
\description{
Given a feature density this function produces two read densities, one for each strand.
}
\usage{
bindDens2readDens(bindDens, fragment, nfrag = 1e+05, bind = 147, 
	minLength = 150, maxLength = 180, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bindDens}{
Numeric vector with the feature density for one chromosome.
}
  \item{fragment}{
Function giving the fragment length distribution.
}
  \item{nfrag}{
Number of fragments that should be simulated to generate the read distribution.
}
  \item{bind}{
Length of binding site.
}
  \item{minLength}{
Minimum fragment length.
}
  \item{maxLength}{
Maximum fragment length.
}
  \item{\dots}{
Further arguments to \code{fragment}.
}
}
\value{
A two column matrix. The first column contains the read density for the forward strand, 
the second column the read density for the reverse strand.
}

\author{
Peter Humburg
}

\seealso{
\code{\link{feat2dens}}, \code{\link{sampleReads}}
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- placeFeatures(start=200, length=1e5)

## calculate feature density
featureDens <- feat2dens(features, length=1e5)

## convert to read density
readDens <- bindDens2readDens(featureDens, fragDens, meanLength=160) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
