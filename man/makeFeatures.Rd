\name{makeFeatures}
\Rdversion{1.1}
\alias{makeFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating a list of genomic features
}
\description{
This function generates a list of genomic features for a single chromosome based on a Markov model.
}
\usage{
makeFeatures(generator = defaultGenerator(), 
	transition = defaultTransition(), init = defaultInit(), 
	start = 1000, length, control = list(), 
	globals = list(minDist = 175), lastFeat = defaultLastFeat(), 
	experimentType = "NucleosomePosition", 
	truncate = FALSE, maxTries = 10, force=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{generator}{
	A named list providing functions to generate the parameters associated with each type of feature. 
	The name of each element in the list is the name of the state the function should be associated with. 
}
  \item{transition}{
	Named list of transition probabilities. Each element is a (named) numeric vector giving the transition
	probabilities for the state indicated by the element's name, i.e., each element of the list is a row
	of the transition probability matrix but zero probabilities can be omitted.
}
  \item{init}{
	Named numeric vector of initial state probabilities. The names have to correspond to state names of the model.
	Zero probabilities may be omitted.
}
  \item{start}{
	Numeric value indicating the position at which the first feature should be placed.
}
  \item{length}{
	Maximum length of DNA that should be covered with features. 
}
  \item{control}{
	Named list with additional arguments to generator functions (one list per generator). 
	Again the names should be the same as the state names.
}
  \item{globals}{
	List of global arguments to be passed to all generator functions.
}
  \item{lastFeat}{
	Named logical vector indicating for each feature type whether it can be the last feature.
}
  \item{experimentType}{
	Type of experiment the simulated features belong to. This is used as the class of the return value.
}
  \item{truncate}{
	Logical value indicating whether the final feature should be truncated to ensure that total length does
	not exceed \code{length} (if \code{FALSE}, a feature that would be truncated is removed instead).
}
  \item{maxTries}{
	Maximum number of attempts made to generate a valid sequence of features. If no valid sequence is 
	generated during the first \code{maxTries} attempts the function will fail either silently 
	(returning an empty sequence) or raise an error, depending on the value of \code{force}.
}
  \item{force}{
	Logical indicating whether the function should be forced to return a feature sequence, even if no
	valid sequence was found. If this is \code{TRUE} an empty sequence will be returned in that case 
	otherwise an error is raised.
}
}
\details{
	This function will generate features from any first order Markov model specified by \code{init}, \code{transition} 
	and \code{generator}. If \code{force} is \code{FALSE} the returned feature sequence is guaranteed to contain at least
	one feature and end in a state that is indicated as possible end state in \code{lastFeat}. Note that the states for
	which \code{lastFeat} is \code{TRUE} are not end states in the sense that the chain is terminated once the state is 
	entered or that the chain remains in the state once it is first entered. Instead this is a mechanism to ensure that 
	some states cannot be the last in the sequence.
	
	Due to the constrains on the total length  of DNA covered by features as well as the possible constraint on the final
	feature of the sequence it is possible to specify models that cannot produce a legal sequence. In other cases it may
	be possible but unlikely to produce a feature sequence that satisfies both constraints. A serious attempt is made to 
	satisfy both requirement, generating a new feature sequence or truncating an existing one if necessary. To ensure that
	this terminates eventually the number of attempts to generate a valid sequence are limited to \code{maxTries}.
	
	In some cases it may be desirable to carry out some post-processing of the feature sequence to ensure that parameters
	of neighbouring features are compatible in some sense. For the default nucleosome positioning simulation the
	function \code{\link{reconcileFeatures}} provides this functionality and \code{\link{placeFeatures}} is an interface
	to \code{makeFeatures} that utilises \code{\link{reconcileFeatures}}. 
}
\value{
	A list of features (with class determined by \code{experimentType}). 
	Each feature is represented by a list of parameters and has a class with the same name as the
	state that generated the feature. In addition all features are of class \code{SimulatedFeature}.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
	Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Functions to generate default values for some of the arguments:
\code{\link{defaultGenerator}}, \code{\link{defaultInit}}, \code{\link{defaultTransition}}, \code{\link{defaultLastFeat}}.

Use \code{\link{feat2dens}} to convert a feature sequence into feature densities.

\code{\link{placeFeatures}} is an interface to \code{makeFeature} for nucleosome positioning.
}
\examples{
set.seed(1)
## generate a (relatively short) sequence of nucleosome features
features <- makeFeatures(length=1e6)

## check the total length of the features
sum(sapply(features, "[[", "length")) ## 995020 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
