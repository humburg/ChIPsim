\name{defaultErrorProb}
\Rdversion{1.1}
\alias{defaultErrorProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replacement probabilities for sequencing errors
}
\description{
For each nucleotide this function provides probabilities indicating how likely it is to be
replaced by any of the other nucleotides should a sequencing error occur.
}
\usage{
defaultErrorProb()
}
%- maybe also 'usage' for other objects documented here.
\details{
The probabilities used here are the ones determined by Dohm \emph{et al.} for \emph{Beta vulgaris}.
They should be more appropriate than a uniform distribution but the use of species specific
rates is recommended where available.
}
\value{
A list of four vectors with replacement probabilities for each nucleotide.
}
\references{
Juliane C. Dohm, Claudio Lottaz, Tatiana Borodina, and Heinz Himmelbauer. Substantial biases in ultra-short
read data sets from high-throughput DNA sequencing. \emph{Nucl. Acids Res.}, pages gkn425+, July 2008.
}
\author{
Peter Humburg
}

\examples{
defaultErrorProb()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
