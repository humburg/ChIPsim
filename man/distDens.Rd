\name{distDens}
\Rdversion{1.1}
\alias{distDens}
\alias{fragDens}
\alias{stableDens}
\alias{indNuc}
\alias{noNuc}
\alias{phaseNuc}
\alias{bindLocDens}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing densities for nucleosome positioning simulation
}
\description{
These functions compute nucleosome densities for a given parameter set (usually provided through one of the feature classes).
}
\usage{
distDens(x, minDist = 175, varDist = 337.5, meanDist = 200)
fragDens(x, minLength, maxLength, meanLength, bind)
indNuc(meanDist = 200, length = 2000, weight = 1)
noNuc(length, weight = 1)
stableDens(x, shift = 10, ratio = 1, weight = 1, stability = 1)
phaseNuc(stable, dist, minDist = 175, length = 2000, meanDist = 200, 
	varDist = (meanDist - minDist) + (meanDist - minDist)^2/2, 
	shift = 10, ratio = 1, weight = 1, stability = 1)
bindLocDens(x, fragLength)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Position at which the density should be evaluated.
}
  \item{minDist}{
Minimum distance between nucleosomes.
}
  \item{varDist}{
Variance of nucleosome distances.
}
  \item{meanDist}{
Mean distance of nucleosomes.
}
  \item{minLength}{
Minimum fragment length.
}
  \item{maxLength}{
Maximum fragment length.
}
  \item{meanLength}{
Mean fragment length.
}
  \item{bind}{
Position of binding site within fragment.
}
  \item{length}{
Length of region.
}
  \item{weight}{
Weight of feature.
}
  \item{stable}{
Density function for stable nucleosome.
}
  \item{dist}{
Density function of distances between nucleosomes.
}
  \item{shift}{
Distance between alternative position for stable nucleosome.
}
  \item{ratio}{
Ratio of probability mass associated with central and alternative positions for stable nucleosome.
}
  \item{stability}{
Stability of stable nucleosome.
}
  \item{fragLength}{
 Length of DNA fragment. If \code{x} is not in [0, 1] this is used to normalize \code{x}.
  }
}

\value{
Density evaluated at the given position. 
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{feat2dens}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{internal}
