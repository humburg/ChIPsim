\name{writeReads}
\Rdversion{1.1}
\alias{writeReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create fastq file from read positions
}
\description{
This is an interface to \code{\link{pos2fastq}} that writes all reads for a given genome to a single file.
}
\usage{
writeReads(readPos, readNames, sequence, quality, file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readPos}{
List of read positions with each component holding the read positions for one chromosome, which are themselves 
two component lists that provide forward and reverse strand positions.
}
  \item{readNames}{
List of the same shape as \code{readPos} providing read names.
}
  \item{sequence}{
Genome reference sequence (a \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}).
}
  \item{quality}{
Read quality scores (see Details).
}
  \item{file}{
Output file.
}
  \item{\dots}{
Further arguments to \code{\link{pos2fastq}}.
}
}
\details{
If \code{quality} looks like it might refer to a fastq file an attempt is made to create a 
\code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object. The read qualities of any 
\code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object passed as \code{quality} 
(directly or indirectly as file name) are extracted and passed on to \code{\link{pos2fastq}} 
as \code{quality} argument. Otherwise it is passed on unchanged. 
}
\value{
The name of the output file.
}
\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pos2fastq}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
