\name{defaultControl}
\Rdversion{1.1}
\alias{defaultControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default parameters for simChIP
}
\description{
Produces a list of parameters for each of the functions used to carry out the various stages of the simulation.
}
\usage{
defaultControl(features = list(), bindDensity = list(), 
	readDensity = list(), readNames = list(), readSequence = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
Parameters for feature generation.
}
  \item{bindDensity}{
Parameters for the conversion of feature sequence into binding site densities. 
}
  \item{readDensity}{
Parameters for the conversion of binding site densities into read densities.
Always provides parameters 
\describe{
\item{\code{fragment}}{Default: \code{\link{fragDens}}}
\item{\code{meanLength}}{Default: 160}
}
}
  \item{readNames}{
Parameters for the generation of read names.
}
  \item{readSequence}{
Parameters for the conversion of read positions into read sequences.
Always provides parameters 
\describe{
\item{\code{qualityFun}}{\code{\link{readQualitySample}}}
\item{\code{errorFun}}{\code{\link{readError}}}
\item{\code{readLen}}{36}
}
}
}
\details{
Any parameters passed as part of list to one of the arguments of \code{defaultControl} will be passed on to
the corresponding function in \code{\link{simChIP}}. The build-in defaults can be overwritten by providing a 
list entry with the same name.
}
\value{
List of parameters for use as the \code{control} argument to \code{\link{simChIP}}.
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{defaultFunctions}}, \code{\link{simChIP}}
}
\examples{
defaultControl()
defaultControl(features=list(maxTail=0), readSequence=list(readLen=50))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

