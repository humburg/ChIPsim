\name{readQualitySample}
\Rdversion{1.1}
\alias{readQualitySample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample read qualities from a list
}
\description{
Given a read sequence and a list of read quality scores this function returns a (possibly truncated)
quality score of the same length as the read.
}
\usage{
readQualitySample(read, qualities, checkLength = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{read}{
A sequence read.
}
  \item{qualities}{
List of sequence read quality scores.
}
  \item{checkLength}{
Flag indicating whether the length of quality scores should be checked to ensure that they are at least as long
as the read. If \code{qualities} contains entries shorter than \code{read} this has to be \code{TRUE}, but see below.  
}
  \item{\dots}{
Further arguments, currently not used.
}
}
\details{
Using \code{checkLength = TRUE} leads to a substantial decrease in performance and is impractical for a large simulation.
To avoid this slow down it is recommended to remove short sequences from \code{qualities} beforehand so that 
\code{checkLength = FALSE} can be used.  
}
\value{
An read quality score string of the same length as \code{read}. 
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{utilities}
