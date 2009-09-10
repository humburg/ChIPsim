\name{FeatureGenerators}
\Rdversion{1.1}
\alias{stableFeature}
\alias{phasedFeature}
\alias{fuzzyFeature}
\alias{nfrFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Features
}
\description{
These functions are used to generate the parameters for different nucleosome positioning related features.
}
\usage{
stableFeature(start, minDist = 175, weight = seq(0.1, 1, by = 0.1), 
	shift = c(0, 5, 10), ratio = seq(0, 4, by = 0.25), 
	stability = seq(0.1, 5, by = 0.1), weightProb, shiftProb, 
	ratioProb, stabilityProb, ...)
phasedFeature(minDist = 175, length = seq(1000, 10000, by = minDist), 
	meanDist = minDist:300, lengthProb, meanDistProb, start, ...)
fuzzyFeature(start, length = seq(1000, 10000, by = 1000), 
	meanDist = 175:400, lengthProb, meanDistProb, ...)
nfrFeature(start, length = seq(50, 500, by = 10), 
	weight = seq(0.1, 1, by = 0.1), lengthProb, weightProb, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
Start location of feature on chromosome.
}
  \item{minDist}{
Minimum distance between nucleosomes.
}
  \item{length}{
A numeric vector giving possible values for the length of the feature.
}
  \item{meanDist}{
A numeric vector giving possible values for the mean distance between nucleosomes.
}
  \item{weight}{
A numeric vector giving possible values for the weight of the feature.
}
  \item{shift}{
A numeric vector giving possible values for the distance between favoured positions of stable nucleosomes.
}
  \item{ratio}{
A numeric vector giving possible values for the ratio between probabilities for 
alternative and central position of stable nucleosomes.
}
  \item{stability}{
A numeric vector giving possible values for the stability of stable nucleosomes.
}
  \item{lengthProb}{
Length distribution of feature. This corresponds to the state duration distribution of the underlying generating model.
}
\item{meanDistProb}{
Distribution of mean distances between nucleosomes.
}
  \item{weightProb}{
Distribution of feature weights.
}
  \item{shiftProb}{
Distribution of distances between favoured positions of stable nucleosome.
}
  \item{ratioProb}{
Ratio distribution.
}
  \item{stabilityProb}{
Stability distribution. 
}
  \item{\dots}{
Further arguments (currently ignored).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of parameters for the corresponding feature type. These parameters are later used to compute nucleosome densities. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simChIP}}
}
\examples{
feature <- stableFeature(200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}

