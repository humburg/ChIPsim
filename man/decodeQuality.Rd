\name{decodeQuality}
\Rdversion{1.1}
\alias{decodeQuality}
\alias{encodeQuality}
\title{
Conversion between numerical and ASCII representation of read qualities
}
\description{
These functions convert an ASCII encoded sequence of read qualities into a numeric vector of error 
probabilities and vice versa.
}
\usage{
decodeQuality(quality, type = c("Illumina", "Sanger", "Solexa"))
encodeQuality(quality, type = c("Illumina", "Sanger", "Solexa"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quality}{
For \code{decodeQuality} a character string representing the read qualities for a single sequence read.
For \code{encodeQuality} a numeric vector of error probabilities.
}
  \item{type}{
Type of encoding to use.
}
}
\details{
See \code{\link{extractQuality}} for a description of the currently supported encodings.
}
\value{
Either a numeric vector of error probabilities or a character string of encoded read quality scores. 
Each entry in the vector corresponds to one character in the input. 
}

\author{
Peter Humburg
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{extractQuality}}
}
\examples{
## decodeQuality and encodeQualty are the inverse operations 
## of each other as one might expect
quality <- "IIIIIIIIIIIIICIIGIIIIGII95III6II-II0"
errorProb <- decodeQuality(quality, type="Sanger")
qualitySanger <- encodeQuality(errorProb, type="Sanger")
all.equal(quality, qualitySanger)

## They can also be used to convert between encodings
qualityIllumina <- encodeQuality(errorProb, type="Illumina")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
