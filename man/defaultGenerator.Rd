\name{defaultGenerator}
\Rdversion{1.1}
\alias{defaultGenerator}
\alias{defaultInit}
\alias{defaultLastFeat}
\alias{defaultTransition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Defaults for Feature Generator
}
\description{
Functions to generate defaults for \code{\link{makeFeatures}}.
}
\usage{
defaultGenerator()
defaultTransition()
defaultInit(prob=c(0.2, 0.05, 0, 0.25, 0.5), 
	states=c("ReversePhasedFeature", "StableFeature", 
		"PhasedFeature", "NFRFeature", "FuzzyFeature"))
defaultLastFeat(isEnd = c(FALSE, rep(TRUE, 4)), 
	states = c("ReversePhasedFeature", "StableFeature", 
		"PhasedFeature", "NFRFeature", "FuzzyFeature"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{prob}{Numeric vector giving the initial state distribution. This will be normalised if the 
		probabilities do not add up to 1.}
	\item{isEnd}{Logical vector indicating which states, i.e. features, are allowed to be last in the sequence.}
	\item{states}{Character vector of state names.}
}

\details{
These functions generate data structures that can be passed as arguments to \code{\link{makeFeatures}}.
Using this set of functions will create a nucleosome positioning simulation. Some of the defaults
can be modified by passing different values to \code{defaultInit} and \code{defaultLastFeat}.
}
\value{
Return values are suitable as arguments \code{generator}, \code{transition}, \code{init} and \code{lastFeat} of 
\code{\link{makeFeatures}}. See the documentation of \code{\link{makeFeatures}} for more detail.
}

\author{
Peter Humburg
}

\examples{
	set.seed(1)
	## generate defaults
	generator <- defaultGenerator()
	transition <- defaultTransition()
	lastFeat <- defaultLastFeat()
	
	## change the initial state distribution such that it 
	## always starts with a fuzzy feature
	init <- defaultInit(c(0, 0, 0, 0, 1))
	
	## now generate some features for a stretch of 1 million base pairs 
	features <- makeFeatures(generator=generator, transition=transition, 
		lastFeat=lastFeat, init=init, length=1e6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{datagen}
