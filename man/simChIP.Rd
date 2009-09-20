\name{simChIP}
\Rdversion{1.1}
\alias{simChIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate ChIP-seq experiments 
}
\description{
This function acts as driver for the simulation. It takes all required arguments and passes
them on to the functions for the different stages of the simulation. The current
defaults will simulate a nucleosome positioning experiment. 
}
\usage{
simChIP(nreads, genome, file, functions = defaultFunctions(), 
	control = defaultControl(), verbose = TRUE, load = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nreads}{
Number of reads to generate.
}
  \item{genome}{
An object of class 'DNAStringSet' or the name of a fasta file containing the genome sequence.
}
  \item{file}{
Base of output file names (see Details).
}
  \item{functions}{
Named list of functions to use for various stages of the simulation, expected names are:
\sQuote{features}, \sQuote{bindDens}, \sQuote{readDens}, \sQuote{sampleReads}, \sQuote{readNames}, \sQuote{readSequence}
}
  \item{control}{
Named list of arguments to be passed to simulation functions (one list per function).
}
  \item{verbose}{
Logical indicating whether progress messages should be printed.
}
  \item{load}{
Logical indicating whether an attempt should be made to load intermediate results from a previous run.
}
}
\details{
The simulation consists of six of stages:
\enumerate{
\item generate feature sequence (for each chromosome): chromosome length -> feature sequence (list)
\item compute binding site density: feature sequence -> binding site density (vector)
\item compute read density: binding site density -> read density (two column matrix, one column for each strand)
\item sample read start sites: read density -> read positions (list)
\item create read names: number of reads -> unique names
\item obtain read sequence and quality: read positions, genome sequence, [qualities] -> output file
}

After each of the first three stages the results of the stage are written to a file and can be reused later. 
File names are created by appending \sQuote{\code{_features.rdata}}, \sQuote{\code{_bindDensity.rdata}} and
\sQuote{\code{_readDensity.rdata}} to \code{file} respectively. Previous results will be loaded for reuse if
\code{load} is \code{TRUE} and files with matching names are found. This is useful to sample repeatedly 
from the same read density or to recover partial results from an interrupted run.

The creation of files can be prevented by setting \code{file = }\dQuote{}. In this case all results will be
returned in a list at the end. Note that this will require more memory since all intermediate results have
to be held until the end. 

The behaviour of the simulation is mainly controlled through the \code{functions} and \code{control} arguments.
They are expected to be lists of the same length with matching names. The names indicate the stage of the simulation
for which the function should be used; elements of \code{control} will be used as arguments for the corresponding 
functions.
}
\value{
A list. The components are typically either lists (with one component per chromosome) or file names
but note that this may depend on the return value of functions listed in \code{functions}.
The components of the returned list are: 
\item{features}{Either a list of generated features or the name of a file containing that list;}
\item{bindDensity}{Either a list with binding site densities or the name of a file containing that list;}
\item{readDensity}{Either a list of read densities or the name of a file containing that list;}
\item{readPosition}{Either a list of read start sites or the name of a file containing that list;}
\item{readSequence}{The return value of the function listed as \sQuote{\code{readSequence}}. 
The default for this the name of the fastq file containing the read sequences;}
\item{readNames}{Either a list of read names or the name of a file containing that list.}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{defaultFunctions}}, \code{\link{defaultControl}}
}
\examples{
\dontrun{
## To run the default nucleosome positioning simulation 
## we can simply run something like the line below.
## This will result in 10 million reads sampled from the genome.
## Of course the file names have to be changed as appropriate. 
simChIP(1e7, genome = "reference.fasta", file = "output/sim_10M")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}

