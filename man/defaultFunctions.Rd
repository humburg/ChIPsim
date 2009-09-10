\name{defaultFunctions}
\Rdversion{1.1}
\alias{defaultFunctions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default functions for simChIP
}
\description{
Provides default functions to carry out the different stages of the ChIP-seq simulation.
}
\usage{
defaultFunctions()
}
%- maybe also 'usage' for other objects documented here.

\value{
A list with components
\item{features}{\code{\link{placeFeatures}}}
\item{bindDensity}{\code{\link{feat2dens}}}
\item{readDensity}{\code{\link{bindDens2readDens}}}
\item{sampleReads}{\code{\link{sampleReads}}}
\item{readSequence}{\code{\link{writeReads}}}
\item{readNames}{\code{\link{simpleNames}}}

}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simChIP}}
}
\examples{
defaultFunctions()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

