\name{writeFASTQ}
\Rdversion{1.1}
\alias{writeFASTQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write read sequences and qualities to a FASTQ formatted file
}
\description{
This is intended to produce the final output of the simulation by providing a fastq file that may 
then be used in an analysis pipeline.
}
\usage{
writeFASTQ(read, quality, name, file, append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read}{
List of read sequences.
}
  \item{quality}{
List of read quality scores.
}
  \item{name}{
Read names.
}
  \item{file}{
File name. If this is \dQuote{} results will be printed to the standard output connection.
}
  \item{append}{
Logical indicating the reads should be appended to an existing file.
}
}
\details{
The first three arguments should have the same length but will be recycled if necessary. 
}
\value{
Called for its side effect.
}

\author{
Peter Humburg 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readSequence}}, \code{\link{readQualitySample}}, \code{\link{writeReads}}
}
\examples{
set.seed(1)

## generate sequence read and quality
quality <- paste(sample(unlist(strsplit(rawToChar(as.raw(33:126)),"")), 
	36, replace = TRUE), collapse="")
read <- paste(sample(c("A", "C", "G", "T"), 36, replace = TRUE), collapse = "")

## write a fastq record
writeFASTQ(read, quality, "read_1", file="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
