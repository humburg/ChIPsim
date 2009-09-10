\name{pos2fastq}
\Rdversion{1.1}
\alias{pos2fastq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert read positions to fastq records
}
\description{
Convert read positions for a single chromosome (both strands) into read sequences + qualities and
write them to file
}
\usage{
pos2fastq(readPos, names, quality, sequence, qualityFun, errorFun, 
	readLen = 36, file, 
	qualityType = c("Illumina", "Sanger", "Solexa"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readPos}{
A list of two numeric vectors (one per strand)
}
  \item{names}{
List of names to use for reads in fastq file. Has to be of same shape as \code{name}.
}
  \item{quality}{
Passed on as argument to \code{qualityFun}.
}
  \item{sequence}{
Reference sequence (a \code{\link[Biostrings]{DNAString}} object). 
}
  \item{qualityFun}{
Function to generate quality scores.
}
  \item{errorFun}{
Function to introduce sequencing errors.
}
  \item{readLen}{
Read length to generate.
}
  \item{file}{
Output file (either file name or connection).
}
  \item{qualityType}{
Encoding to use for read quality scores.
}
  \item{\dots}{
Further arguments (see Details).
}
}
\details{
Arguments passed as part of \dots will be passed on to \code{qualityFun}, except an argument called \code{prob} which is
passed on to errorFun instead if present.
}
\value{
Invisibly returns the number of records that were written.
}
\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{readError}} for a possible choice of \code{errorFun} and \code{\link{readQualitySample}} for a simple 
\code{qualityFun}.
}
\examples{
set.seed(1)

## a function to generate random read qualities (in Sanger format)
randomQuality <- function(read, ...){
	paste(sample(unlist(strsplit(rawToChar(as.raw(33:126)),"")), 
		length(read), replace = TRUE), collapse="")
} 

## generate a reference sequence
chromosome <- DNAString(paste(sample(c("A", "C", "G", "T"), 
	1e5, replace = TRUE), collapse = ""))

## and a few read positions
reads <- list(sample(100:9900, 5), sample(100:9900, 5))
names <- list(paste("read", 1:5, sep="_"), paste("read", 6:10, sep="_"))

## convert to fastq format
pos2fastq(reads, names, sequence=chromosome, qualityFun=randomQuality, 
	errorFun=readError, file="")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
