\name{extractQuality}
\Rdversion{1.1}
\alias{extractQuality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain read qualities from a Fastq file or ShortReadQ object
}
\description{
Converts the read qualities encoded in fastq formatted files into error probabilities.
}
\usage{
extractQuality(reads, minLength = 25, dir, 
	type = c("Illumina", "Sanger", "Solexa"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{
Either the name of a fastq file or a \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object (see Details).
}
  \item{minLength}{
Minimum read length required.
}
  \item{dir}{
Directory of fastq file.
}
  \item{type}{
Character string indicating the format the qualities are encoded in (see Details).
}
}
\details{
If \code{reads} and \code{dir} are character strings it is assumed that \file{dir/reads} is the name 
of a fastq file. Otherwise \code{reads} should be a \code{\link[ShortRead]{ShortReadQ}} object in which
case \code{dir} is ignored.

Currently three different encodings of read qualities are supported. The encoding has to be selected via the 
\code{type} argument. The supported formats are
\describe{
	\item{Illumina}{The format currently used by Illumina (version 1.3). This is a phred score between 0 and 40 
	encoded as ASCII characters 64 to 104. [default]}
	\item{Sanger}{The Sanger format uses a phred quality score between 0 and 93 encoded as ASCII characters 33 to 126.}
	\item{Solexa}{The old Solexa format previously used by Solexa/Illumina uses a quality score between -5 and 40 
	encoded as ASCII characters 59 to 104.}
} 
}
\value{
A list with a vector of error probabilities for each read in \code{reads} that is at least \code{minLength} 
nucleotides long.
}

\author{
Peter Humburg
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{decodeQuality}}, \code{\link{readQualitySample}}
}
\examples{
\dontrun{
	## load reads from a fastq file with Sanger encoding
	qualities <- extractQuality("test.fastq", dir=".", type="Sanger")
	
	## extract error probabilities for first 25bp of each read
	qualities25 <- sapply(qualities, "[", 1:25)
	
	## plot average quality for each position
	plot(rowMeans(qualities25), type='b', xlab="Read position", 
		ylab="Error probability")
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
