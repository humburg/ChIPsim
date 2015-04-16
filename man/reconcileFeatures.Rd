\name{reconcileFeatures}
\Rdversion{1.1}
\alias{reconcileFeatures}
\alias{reconcileFeatures.NucleosomePosition}
\alias{reconcileFeatures.SimulatedExperiment}
\alias{reconcileFeatures.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post-processing of simulated features
}
\description{
The \code{reconcileFeatures} functions provide a facility to post-process a list of features representing
a simulated experiment. \code{reconcileFeatures} is an S3 generic, new functions can be added for additional 
types of experiment. The current default is to call \code{reconcileFeatures.SimulatedExperiment} which, if
called without further arguments, will simply return the feature list unchanged. 
}
\usage{
reconcileFeatures(features, ...)
\method{reconcileFeatures}{default}(features, ...)
\method{reconcileFeatures}{SimulatedExperiment}(features, defaultValues=list(), ...)
\method{reconcileFeatures}{NucleosomePosition}(features, defaultMeanDist = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
List of simulated features.
}
  \item{defaultValues}{
Named list of default parameter values. The method for class \code{SimulatedExperiment} ensures that
all features have at least the parameters listed in \code{defaultValues}, adding them where necessary.
}
  \item{defaultMeanDist}{
Default value for the average distance between nucleosomes for nucleosome positioning experiments.
}
\item{\dots}{Further arguments to future functions.}
}

\value{
A list of features of the same class as \code{features}.
}

\author{
Peter Humburg
}

\seealso{
\code{\link{makeFeatures}}, \code{\link{placeFeatures}}
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- makeFeatures(length=1e6, )

## check the total length of the features
sum(sapply(features, "[[", "length")) ## 995020

## reconcile features to ensure smooth transitions
## For experiments of class NucleosomePosition this 
## also combines some features and introduces
## some overlap between them.
features <- reconcileFeatures(features)
 
## check the total length of the features again
sum(sapply(features, "[[", "length")) ## 984170
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
