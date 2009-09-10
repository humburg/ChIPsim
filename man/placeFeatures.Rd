\name{placeFeatures}
\Rdversion{1.1}
\alias{placeFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating and reconciling a feature sequence
}
\description{
This function provides an interface to \code{\link{makeFeatures}} and \code{\link{reconcileFeatures}} that combines
both steps of the feature generation process.
}
\usage{
placeFeatures(..., maxTail = 0.01, 
	compoundFeatures=list("StablePhasedFeature"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Arguments to \code{\link{makeFeatures}}.
}
  \item{maxTail}{
Maximum portion of total length of chromosome that may be left free of features (see Details).
}
  \item{compoundFeatures}{
List of feature classes that are produced by combining two features. This may happen during the call 
to \code{\link{reconcileFeatures}} and requires special handling when extending the feature list.  
}
}
\details{
This function (as well as \code{\link{makeFeatures}} which it calls) tries to fill as much of the genomic region
with features as possible, i.e. an attempt is made to produce a feature sequence that covers \code{length}
base pairs. In most cases the sequence will be slightly shorter. The \code{maxTail} argument determines how long
a region without any features at the end of the genomic region is acceptable (as fraction of the total length).
Note however that even \code{maxTail = 0} does not guarantee a feature sequence of exactly the requested length.
}
\value{
A list of simulated features. The class of the return value as well as the features generated depend on
the arguments passed to \code{\link{makeFeatures}}. 
}

\author{
Peter Humburg
}
\note{
Using the \code{\link{reconcileFeatures}} mechanism it is possible to introduce dependence between neighbouring features 
that is not easily expressed in terms of a simple Markov model. In some cases the same effect can be achieved by 
introducing additional states into the model but it may be more convenient to simply post-process the feature sequence.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeFeatures}}, \code{\link{reconcileFeatures}}
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- placeFeatures(length=1e6, maxTail = 0)

## check the total length of the features
sum(sapply(features, "[[", "length")) ## 990509
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
