\name{readError}
\Rdversion{1.1}
\alias{readError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Introduce errors into read sequence based on quality scores
}
\description{
Given a read sequence and quality this function introduces errors by first choosing positions 
that should be modified based on the quality score and then exchanges nucleotides based on the
probabilities given in \code{prob}. 
}
\usage{
readError(read, qual, alphabet = c("A", "C", "G", "T"), 
	prob = defaultErrorProb(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read}{
A character string representing a read sequence.
}
  \item{qual}{
Numeric vector of read error probabilities.
}
  \item{alphabet}{
Alphabet used for read sequence.
}
  \item{prob}{
  Nucleotide exchange probabilities.
}
  \item{\dots}{
Further arguments (currently ignored).
}
}
\details{
If the read sequence contains letters that are not part of \code{alphabet} they are replaced 
by the first entry of \code{alphabet} before positions of sequencing errors are determined.
The alphabet used has to match the names used in \code{prob}.
}
\value{
The modified read sequence.
}

\author{
Peter Humburg
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{defaultErrorProb}}, \code{\link{readSequence}}
}
\examples{
set.seed(42)

## generate sequence read and quality
quality <- paste(sample(unlist(strsplit(rawToChar(as.raw(33:126)),"")), 
	36, replace = TRUE), collapse="")
errorProb <- decodeQuality(quality, type = "Sanger")
read <- paste(sample(c("A", "C", "G", "T"), 36, replace = TRUE), 
	collapse = "")

## use readError to introduce sequencing errors
read2 <- readError(read, errorProb)

all.equal(read, read2)  ## "1 string mismatch"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
