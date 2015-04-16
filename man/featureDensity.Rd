\name{featureDensity}
\Rdversion{1.1}
\alias{featureDensity}
\alias{featureDensity.FuzzyFeature}
\alias{featureDensity.NFRFeature}
\alias{featureDensity.ReversePhasedFeature}
\alias{featureDensity.StableFeature}
\alias{featureDensity.StablePhasedFeature}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing density for a given feature
}
\description{
This set of functions is used to generate the density of individual features of different types. 
\code{featureDensity} is an S3 generic, functions may be defined for different feature classes. 
}
\usage{
featureDensity(x, ...)
\method{featureDensity}{StableFeature}(x, stable=stableDens, background=FALSE, ...)
\method{featureDensity}{StablePhasedFeature}(x, stable=stableDens, dist=distDens, background=FALSE, ...)
\method{featureDensity}{ReversePhasedFeature}(x, stable=stableDens, dist=distDens, background=FALSE, ...)
\method{featureDensity}{NFRFeature}(x, background=FALSE, ...)
\method{featureDensity}{FuzzyFeature}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The feature for which the density should be computed.
}
  \item{stable}{
Function that should be used to compute the density of a stable feature.
}
  \item{dist}{
Function that should be used to compute the distribution of distances between adjacent features.
}
  \item{background}{
Logical indicating whether uniform background should be added to the feature.
}
  \item{\dots}{
Arguments to future functions.
}
}
\details{
 These functions are used internally by \code{\link{feat2dens}}. There should be no need to call them 
 directly but it is important to supply suitable \code{featureDensity} functions for new feature types. 
}
\value{
A two column matrix. The first column contains the density, the second the weight at each position.
}
\author{
Peter Humburg
}
\seealso{
\code{\link{feat2dens}}, \code{\link{makeFeatures}}
}
\examples{
## Create a single StableFeature
feature <- stableFeature(start = 200, weight = 0.8, shift = 10, 
	stability = 1, ratio = 1)

## Convert the feature into a density (without background)
featDens <- featureDensity(feature)

## If we want featureDensity to automatically add uniform background
## we have to ensure that the feature has a 'meanDist' parameter
## (this is usually added by 'reconcileFeatures').
feature$meanDist <- 200
featDens2 <- featureDensity(feature, background = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
