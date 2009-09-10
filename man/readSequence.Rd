\name{readSequence}
\Rdversion{1.1}
\alias{readSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert read position into read sequence
}
\description{
Given a read position, a reference sequence, a strand and a read length this function returns 
the read sequence.
}
\usage{
readSequence(readPos, sequence, strand, readLen = 36)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readPos}{
Numeric value indicating the start position on the chromosome.
}
  \item{sequence}{
Chromosome sequence (a \code{\link[Biostrings]{DNAString}})
}
  \item{strand}{
Strand indicator (+1 / -1)
}
  \item{readLen}{
Length of read.
}
}

\value{
Read sequence.
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readError}}, \code{\link{writeReads}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
