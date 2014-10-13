\name{simpleNames}
\Rdversion{1.1}
\alias{simpleNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate unique read names
}
\description{
Generates a set of unique (and very simple) read names.
}
\usage{
simpleNames(n, nameBase = "read")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of names to generate.
}
  \item{nameBase}{
Base name to use.
}
}
\value{
A character vector with entries of the form \sQuote{nameBase_i} where i runs from 1 to \code{n}.
}
\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
simpleNames(5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
