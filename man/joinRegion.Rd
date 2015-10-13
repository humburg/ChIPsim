\name{joinRegion}
\Rdversion{1.1}
\alias{joinRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combining two feature densities
}
\description{
Function to take two vectors of feature densities and combine them into a single
vector, using overlap between the two densities and smoothing the transition. 
}
\usage{
joinRegion(left, right, overlap, overlapWeights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{left}{
First density vector.
}
  \item{right}{
Second density vector.
}
  \item{overlap}{
Overlap between the two features.
}
  \item{overlapWeights}{
Weights for overlapping region.
}
}

\value{
Returns the combined density vector.
}

\author{
Peter Humburg
}
\note{
This function is used as part of \code{\link{feat2dens}} and there should be no need to call it directly although it
may be useful for possible extensions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
