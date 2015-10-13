\name{feat2dens}
\Rdversion{1.1}
\alias{feat2dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a list of features into a feature density 
}
\description{
Given a list of features (as produced by \code{\link{makeFeatures}}) computes the feature density for each
and combines them into a chromosome wide density.
}
\usage{
feat2dens(features, length, featureBgr = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
A list of features. 
}
  \item{length}{
Total length of feature density vector (i.e. chromosome length). If this is missing the length is inferred from the feature
parameters. 
}
  \item{featureBgr}{
Logical indicating whether feature specific background should be added to the density. If this is \code{TRUE} 
the resulting density for each feature is a mixture of the feature density and a fuzzy, i.e. uniform, 
feature density. The weights of the components are determined by the feature weight.
}
  \item{\dots}{
Further arguments to \code{\link{featureDensity}}.
}
}

\value{
A vector with the feature density for each position along the chromosome.
}

\author{
Peter Humburg
}

\seealso{
The majority of the work is done by calls to \code{\link{featureDensity}} and \code{\link{joinRegion}}.
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- placeFeatures(start=200, length=1e5)

## calculate density
featureDens <- feat2dens(features, length=1e5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
