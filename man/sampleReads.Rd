\name{sampleReads}
\Rdversion{1.1}
\alias{sampleReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling sequence read positions from a read density.
}
\description{
Given a read density this function returns the starting positions of sequence reads. 
}
\usage{
sampleReads(readDens, nreads = 6e+06, strandProb = c(0.5, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readDens}{
A two column matrix of read densities (as produced by \code{\link{bindDens2readDens}}).
}
  \item{nreads}{
Number of read positions to generate.
}
  \item{strandProb}{
A numeric vector with two elements giving weights for forward and reverse strand.  
}
}
\details{
The expected number of reads for each strand is \code{strandProb * nreads}.
}
\value{
A list with components \code{fwd} and \code{rev} giving the read positions on the 
forward and reverse strand respectively.
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bindDens2readDens}}
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- placeFeatures(start=200, length=1e5)

## calculate feature density
featureDens <- feat2dens(features, length=1e5)

## convert to read density
readDens <- bindDens2readDens(featureDens, fragDens, meanLength=160) 

## sample reads
## of course you would usually want a much larger number
readPos <- sampleReads(readDens, nreads=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
