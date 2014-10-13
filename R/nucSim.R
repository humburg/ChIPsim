## functions for the simulation of nucleosome ChIP-seq

############################# SECTION: generic functions ######################################################
featureDensity <- function(x, ...){
	UseMethod("featureDensity")
}

reconcileFeatures <- function(features, ...)
	UseMethod("reconcileFeatures")

############################# SECTION: Defaults ###############################################################

reconcileFeatures.SimulatedExperiment <- function(features, defaultValues = list(), ...){
	features <- lapply(features, function(f){
				if(length(defaultValues) > 0)
					for(param in names(defaultValues)){
						if(is.null(f[[param]])) f[[param]] <- defaultValues[[param]]
					}
				f$overlap <- 0
				currentClass <- class(f)
				class(f) <- c(currentClass[-length(currentClass)], 
						"ReconciledFeature", currentClass[length(currentClass)])
				f
			})
	features
}

reconcileFeatures.default <- function(features, ...) reconcileFeatures.SimulatedExperiment(features, ...)

############################# SECTION: generating feature sequence ############################################

## default parameters for feature generator
defaultGenerator <- function(){
	list(StableFeature=stableFeature, PhasedFeature=phasedFeature, 
			ReversePhasedFeature=phasedFeature, NFRFeature=nfrFeature, FuzzyFeature=fuzzyFeature)
}

defaultTransition <- function(){
	trans <- list(StableFeature=c(0.6, 0.4), PhasedFeature=c(0.2,  0.3, 0.5), 
			ReversePhasedFeature=1, NFRFeature=c(0.5, 0.3, 0.2), FuzzyFeature=c(0.5, 0.5))
	names(trans$StableFeature) <- c("PhasedFeature", "NFRFeature")
	names(trans$PhasedFeature) <- c("ReversePhasedFeature", "NFRFeature", "FuzzyFeature")
	names(trans$ReversePhasedFeature) <- "StableFeature"
	names(trans$NFRFeature) <- c("StableFeature", "ReversePhasedFeature", "FuzzyFeature")
	names(trans$FuzzyFeature) <- c("ReversePhasedFeature", "NFRFeature")
	
	trans <- lapply(trans, "class<-", "StateDistribution")
	trans
}

defaultInit <- function(prob=c(0.2, 0.05, 0, 0.25, 0.5), 
		states=c("ReversePhasedFeature", "StableFeature", "PhasedFeature", "NFRFeature", "FuzzyFeature")){
	if(sum(prob) != 1) prob <- prob/sum(prob)
	names(prob) <- states 
	class(prob) <- "StateDistribution"
	prob
}

defaultLastFeat <- function(isEnd = c(FALSE, rep(TRUE, 4)), 
		states = c("ReversePhasedFeature", "StableFeature", "PhasedFeature", "NFRFeature", "FuzzyFeature")){ 
	if(!any(isEnd))
		stop("No end states specified.")
	if(!is.null(states))
		names(isEnd) <- states
	if(is.null(names(isEnd))) 
		stop("State names are required.")
	isEnd
}
	

## Markov chain based feature generation
## generator: named list of functions to generate different types of features
## transition: named list of transition probabilities, each entry is expected to be a named numeric vector
##             corresponding to the non zero entries in a row of the transition matrix
## init: vector with initial probabilities for each feature type
## start: start position of first feature
## length: total length of genome (or chromosome) to fill with features
## control: named list with additional arguments to generator functions (one list per generator)
## globals: list of global arguments to be passed to all generator functions
## lastFeat: logical vector indicating for each feature type whether it can be the last feature
## truncate: logical indicating whether the final feature should be truncated to ensure that total length does
##			not exceed 'length' (if FALSE, a feature that would would be truncated is removed instead).
## maxTries: maximum number of attempts made to generate a valid sequence of features. If no valid sequence is 
##			generated during the first 'maxTries' attempts the function will fail either silently 
##			(returning an empty sequence) or raise an error, depending on the value of 'force'
## force: logical indicating whether the function should be forced to return a feature sequence, even if no
##			valid sequence was found. If this is TRUE an empty sequence will be returned in that case 
##			otherwise an error is raised.
makeFeatures <- function(generator=defaultGenerator(), transition=defaultTransition(), init=defaultInit(), 
		start=1000, length, control=list(), globals=list(minDist=175), lastFeat=defaultLastFeat(), 
		experimentType = "NucleosomePosition", truncate=FALSE, maxTries=10, force=FALSE){
	
	test <- sapply(generator, function(x) !is.function(x) && !is.name(x))
	if(any(test)) 
		generator[test] <- lapply(generator[test], as.name)

	features <- list()
	if(inherits(init, "StateDistribution")){
		state <- sample(names(init), 1, prob=init)
		curLength <- start
	}
	else if(inherits(init, "SimulatedFeature")){
		features[[1]] <- init
		state <- sample(names(transition[[class(init)[1]]]), 1, prob=transition[[class(init)[1]]])
		curLength <- start + init$length
	}
	else stop("Argument 'init' is of class ", class(init)[1], ". Expected 'StateDistribution' or 'SimulatedFeature'.")
	
	while(curLength < length){
		## get parameters for current state
		pars <- c(list(start=curLength+1), control[[state]])
		pars[names(globals)] <- globals
		## construct and evaluate call to generator function
		feat <- do.call(generator[[state]], pars)
		class(feat) <- c(state, "SimulatedFeature")
		
		curLength <- curLength + feat$length
		if(curLength > length) {
			if(lastFeat[state] && truncate)
				feat$length <- feat$length - curLength + length
			else break
		}
		features[[length(features)+1]] <- feat
		state <- sample(names(transition[[state]]), 1, prob=transition[[state]])
	}
	
	tries <- 1
	## if no features were generated, try again
	while(tries < maxTries && length(features) == 0){
		features <- Recall(generator, transition, init, start, length, control, globals, lastFeat, truncate, 1, TRUE)
		tries <- tries + 1
	}
	if(length(features) == 0){
		if(force) return(list())
		
		helpText <- ""
		if(!truncate)
			helpText <- paste(" This could be because the initial feature was longer than the maximum length (", 
					length - start, "). Consider using 'truncate=TRUE' or adjust your model parameters to ensure",
					" that at least one feature is generated.", sep="")
		else helpText <- " Check your model parameters."
		stop("Failed to generate any features.", helpText)
	}
	## ensure that last feature is valid end state
	validEnd <- names(lastFeat)[lastFeat]
	while(tries < maxTries && !any(sapply(validEnd, inherits, x=features[[length(features)]]))){
		previous <- features[[length(features) - 1]]
		curLength <- curLength - previous$length
		## remove (invalid) last feature and restart MC in what is now the last state
		features <- features[-c(length(features)-1, length(features))]
		features <- c(features, Recall(generator, transition, previous, curLength + 1, 
						length, control, globals, lastFeat, truncate, 1, TRUE))
		tries <- tries + 1		
	}
	## if we still didn't find a valid end state we try to truncate the feature sequence to get valid sequence
	if(!any(sapply(validEnd, inherits, x=features[[length(features)]]))){
		endPositions <- sapply(features, function(x) any(sapply(validEnd, inherits, x=x)))
		if(any(endPositions)){
			features <- features[1:max(which(endPositions))]
			warning("Feature sequence was truncated to ensure that it ends in a valid feature. ",
					"The feature sequence may be much shorter than requested.")
		}
		else{
			if(force) return(list())
			stop("Failed to generate valid end state. Consider adjusting your model parameters.")
		}
	}
	
	class(features) <- c(experimentType, "SimulatedExperiment")
	features	
}


############################# SECTION: placing features #######################################################
## Although all features have a start position this is only used for the first feature. All remaining
## features are placed based on the position, length and overlap of the previous feature

## generating parameters for a stable nucleosome
stableFeature <- function(start, minDist=175, weight=seq(0.1, 1, by=0.1), shift=c(0, 5, 10), 
		ratio=seq(0, 4, by=0.25), stability=seq(0.1, 5, by=0.1), 
		weightProb, shiftProb, ratioProb, stabilityProb, ...){
	if(missing(weightProb)) weightProb <- dbeta(weight, 2, 2)
	if(missing(shiftProb)) shiftProb <- rep(1/length(shift), length(shift))
	if(missing(ratioProb)) ratioProb <- rep(1/length(ratio), length(ratio))
	if(missing(stabilityProb)) stabilityProb <- dgamma(stability, 2, 2)
	
	## choose parameters
	params <- list()
	## general feature parameters
	params$start <- start - minDist                    ## start position (adjusted to include left flank)
	params$length <- 2*minDist + 1                     ## length of window (just the stable nucleosome)
	params$weight <- ifelse(length(weight) > 1, sample(weight, 1, prob=weightProb), weight) 
													   ## weight of feature
	## stable nucleosome specific parameters
	params$minDist <- minDist						   ## this is constant and only recorded here for convenience, could be dropped
	params$shift <- ifelse(length(shift) > 1, sample(shift, 1, prob=shiftProb), shift)  
													   ## distance between central and alternative nucleosome positions
	params$ratio <- ifelse(length(ratio) > 1, sample(ratio, 1, prob=ratioProb), ratio)   
													   ## ratio between alternative and central position probabilities
	params$stability <- ifelse(length(stability) > 1, sample(stability, 1, prob=stabilityProb), stability)
                                                       ## stability of nucleosome
													   
	class(params) <- c("StableFeature", "SimulatedFeature")
	
	params
}

## generating parameters for phased nucleosomes
## Note that this does not include parameters related to the stable nucleosome that induces the phasing
phasedFeature <- function(minDist=175, length=seq(1000, 10000, by=minDist), meanDist=minDist:300, 
		lengthProb, meanDistProb, start, ...){
	if(missing(lengthProb)) lengthProb <- rep(1/length(length), length(length))
	if(missing(meanDistProb)) meanDistProb <- dgamma(meanDist - minDist, 5, scale=15) #rep(1/length(meanDist), length(meanDist))
	
	## choose parameters
	params <- list()
	## general feature parameters
	params$start <- start                                 
	params$length <- ifelse(length(length) > 1, sample(length, 1, prob=lengthProb), length) 
													   ## length of window
	params$weight <- NA                                ## weight of feature is determined by stable nucleosome
	## phasing specific parameters
	params$minDist <- minDist						   ## this is constant and only recorded here for convenience, could be dropped
	params$meanDist <- ifelse(length(meanDist) > 1, sample(meanDist, 1, prob=meanDistProb), meanDist)
												       ## average distance between nucleosomes
	
	class(params) <- c("PhasedFeature", "SimulatedFeature")
	
	params
}

## generating parameters for nucleosome free region
nfrFeature <- function(start, length=seq(50, 500, by=10), weight=seq(0.1, 1, by=0.1), 
		lengthProb, weightProb, ...){
	if(missing(lengthProb)) lengthProb <- rep(1/length(length), length(length))
	if(missing(weightProb)) weightProb <- rep(1/length(weight), length(weight))
	
	## choose parameters
	params <- list()
	## general feature parameters
	params$start <- start                             ## start position
	params$length <- ifelse(length(length) > 1, sample(length, 1, prob=lengthProb), length) 
														## length of window
	params$weight <- ifelse(length(weight) > 1, sample(weight, 1, prob=weightProb), weight) 
														## weight of feature
	
	class(params) <- c("NFRFeature", "SimulatedFeature")
	
	params
}

## generating parameters for independently placed (fuzzy) nucleosomes
fuzzyFeature <- function(start, length=seq(1000, 1e4, by=1000), meanDist=175:400,
		lengthProb, meanDistProb, ...){
	if(missing(lengthProb)) lengthProb <- rep(1/length(length), length(length))
	if(missing(meanDistProb)) meanDistProb <- dgamma(meanDist-min(meanDist), 10, scale=10)
	
	## choose parameters
	params <- list()
	## general feature parameters
	params$start <- start                             ## start position
	params$length <- ifelse(length(length) > 1, sample(length, 1, prob=lengthProb), length) 
													  ## length of window
	params$weight <- 1                                ## weight of feature (fuzzy features are only generated in the absence of other features)
	## fuzzy nucleosome specific parameters
	params$meanDist <- ifelse(length(meanDist) > 1, sample(meanDist, 1, prob=meanDistProb), meanDist)
													  ## average distance between nucleosomes

	class(params) <- c("FuzzyFeature", "SimulatedFeature")
	
	params
}

## given a list of features, ensure parameters of neighbouring features are compatible and determine required
## overlap between features
reconcileFeatures.NucleosomePosition <- function(features, defaultMeanDist=200, ...){
	exClass <- class(features)
	lastMeanDist <- sapply(features, function(x) if(!is.null(x$meanDist)) x$meanDist else NA)
	lastMeanDist <- lastMeanDist[!is.na(lastMeanDist)]
	if(length(lastMeanDist) == 0) lastMeanDist <- defaultMeanDist
	else lastMeanDist <- lastMeanDist[[1]]
	
	i <- 1
	while(TRUE){
		if(inherits(features[[i]], "StableFeature") && !inherits(features[[i]], "ReconciledFeature")){
			switch(class(features[i+1][[1]])[1],
					PhasedFeature={
						## combine stable nucleosome and following phased nucleosomes into single feature
						features[[i]]$length <- features[[i]]$length + features[[i+1]]$length - features[[i+1]]$minDist
						features[[i]]$minDist <- features[[i+1]]$minDist
						features[[i]]$meanDist <- features[[i+1]]$meanDist
						
						features[[i]]$overlap <- NA
						features[[i]]$overlapWeight <- NA
						
						class(features[[i]]) <- unique(c("StablePhasedFeature", 
										class(features[[i]])[-length(class(features[[i]]))],
										class(features[[i+1]])[-length(class(features[[i+1]]))], 
										"ReconciledFeature", "SimulatedFeature"))
						features <- features[-(i+1)]
					},
					NFRFeature={
						## ensure starting position is consistent and determine overlap
						features[[i]]$overlap <- min(features[[i]]$minDist, floor(0.75*features[[i+1]]$length))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						features[[i]]$meanDist <- lastMeanDist
						
						class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))],
								"ReconciledFeature", "SimulatedFeature")
						i <- i + 1
					},
					ifelse(is.null(features[i+1][[1]]), {features[[i]]$overlap <- NA ;break}, 
							stop("Illegal feature transition: ", class(features[[i]])[1], " -> ", class(features[[i+1]])[1]))
			)		
		}
		if(inherits(features[[i]], "PhasedFeature") || inherits(features[[i]], "ReversePhasedFeature")){
			lastMeanDist <- features[[i]]$meanDist
			switch(class(features[i+1][[1]])[1],
					StableFeature={
						## copy parameter values for following stable nucleosome
						features[[i]]$weight <- features[[i+1]]$weight
						features[[i]]$shift <- features[[i+1]]$shift
						features[[i]]$ratio <- features[[i+1]]$ratio
						features[[i]]$stability <- features[[i+1]]$stability
						
						features[[i]]$overlap <- features[[i]]$minDist
						features[[i]]$overlapWeight <- sum
#						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						features[[i+1]]$meanDist <- features[[i]]$meanDist
						
						class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))],
								"ReconciledFeature", "SimulatedFeature")
						i <- i+1
					},
					## since the next feature is not a stable nucleosome this PhasedFeature 
					## really is a StablePhasedFeature already. Just need to determine overlap
					ReversePhasedFeature={
						features[[i]]$overlap <- floor(min(features[[i]]$length*0.25, features[[i+1]]$length*0.25, 1000))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						i <- i+1
					},
					NFRFeature={
						features[[i]]$overlap <- floor(min(features[[i]]$length*0.25, features[[i+1]]$length*0.75, 
								features[[i]]$meanDist))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						i <- i+1
					},
					FuzzyFeature={
						features[[i]]$overlap <- floor(min(features[[i]]$length*0.25, features[[i+1]]$length*0.25, 1000))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						i <- i+1
					},
					ifelse(is.null(features[i+1][[1]]), {features[[i]]$overlap <- NA ;break}, 
							stop("Illegal feature transition: ", class(features[[i]])[1], " -> ", class(features[[i+1]])[1]))
			)
		}
		if(inherits(features[[i]], "NFRFeature")){
			switch(class(features[i+1][[1]])[1],
				StableFeature={
					features[[i]]$overlap <- min(features[[i+1]]$minDist, floor(0.75 * features[[i]]$length))
					features[[i]]$overlapWeight <- sum
					features[[i]]$meanDist <- lastMeanDist
					features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
					
					class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))],
							"ReconciledFeature", "SimulatedFeature")
					i <- i + 1
				},
				FuzzyFeature=,
				ReversePhasedFeature={
					features[[i]]$overlap <- floor(min(features[[i+1]]$length*0.25, features[[i]]$length*0.75, 
							features[[i+1]]$meanDist))
					features[[i]]$overlapWeight <- sum
					features[[i]]$meanDist <- features[[i+1]]$meanDist
					features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
					
					lastMeanDist <- features[[i]]$meanDist
					class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))], "ReconciledFeature",
							"SimulatedFeature")
					i <- i + 1
				},
				ifelse(is.null(features[i+1][[1]]), {features[[i]]$overlap <- NA ;break}, 
						stop("Illegal feature transition: ", class(features[[i]])[1], " -> ", class(features[[i+1]])[1]))
			)
		}
		if(inherits(features[[i]], "FuzzyFeature")){
			lastMeanDist <- features[[i]]$meanDist
			switch(class(features[i+1][[1]])[1],
					ReversePhasedFeature={
						features[[i]]$overlap <- floor(min(features[[i]]$length*0.25, features[[i+1]]$length*0.25, 1000))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))], "ReconciledFeature",
								"SimulatedFeature")
						i <- i + 1
					},
					NFRFeature={
						features[[i]]$overlap <- floor(min(features[[i]]$length*0.25, features[[i+1]]$length*0.75, 
								features[[i]]$meanDist))
						features[[i]]$overlapWeight <- sum
						features[[i+1]]$start <- features[[i]]$start + features[[i]]$length - features[[i]]$overlap
						
						class(features[[i]]) <- c(class(features[[i]])[-length(class(features[[i]]))], "ReconciledFeature",
								"SimulatedFeature")
						i <- i + 1
					},
					ifelse(is.null(features[i+1][[1]]), {features[[i]]$overlap <- NA ;break}, 
							stop("Illegal feature transition: ", class(features[[i]])[1], " -> ", class(features[[i+1]])[1]))
			)
		}
		
	}
	class(features) <- c(exClass[-length(exClass)], "ReconciledSimulatedExperiment", "SimulatedExperiment")	
	features
}

## function to generate and reconcile feature sequence
## all arguments are passed on to makeFeatures
placeFeatures <- function(..., maxTail=0.01, compoundFeatures=list("StablePhasedFeature")){
	features <- makeFeatures(...)
	lastState <- features[[length(features)]]
	features <- reconcileFeatures(features)
	
	## ensure gap at end is not too large
	args <- list(...)
	length <- args$length
	gap <- length - features[[length(features)]]$start - features[[length(features)]]$length
	
	while(gap > maxTail * length){
		args$start <- features[[length(features)]]$start
		if(is.null(args$transition)) args$transition <- formals(makeFeatures)$transition
		args$init <- lastState
		
		gapFeat <- do.call(makeFeatures, args)
		if(any(sapply(compoundFeatures, inherits, x=features[[length(features)]])))
			gapFeat <- reconcileFeatures(c(features[length(features)], gapFeat[-1]))
		else gapFeat <- reconcileFeatures(gapFeat)
		features <- c(features[-length(features)], gapFeat) 
			
		lastGap <- gap
		gap <- length - features[[length(features)]]$start - features[[length(features)]]$length
		if(gap == lastGap) break
	}
	
	if(gap < 0) features[[length(features)]]$length <- features[[length(features)]]$length + gap
	
	features
}


######################### SECTION: Turn features into densities ###############################################
## x: feature object
featureDensity.StablePhasedFeature <- function(x, stable=stableDens, dist=distDens, background=FALSE, ...){
	dens <- phaseNuc(stable, dist, minDist=x$minDist, length=x$length, meanDist=x$meanDist, shift=x$shift, ratio=x$ratio, 
				weight=x$weight, stability=x$stability)
	if(background && x$weight < 1){
		bg <- indNuc(meanDist=x$meanDist, weight=1-x$weight, length=x$length)
		dens[, 1] <- dens[,1] + bg[,1]
	}
	dens
}

featureDensity.ReversePhasedFeature <- function(x, stable=stableDens, dist=distDens, background=FALSE, ...){
	dens <- phaseNuc(stable, dist, minDist=x$minDist, length=x$length+x$minDist, meanDist=x$meanDist, shift=x$shift, 
			ratio=x$ratio, weight=x$weight, stability=x$stability)
	densStable <- stable(-x$minDist:x$minDist, shift=x$shift, ratio=x$ratio, weight=x$weight, stability=x$stability)
	dens[1:(2*x$minDist+1), 1] <- dens[1:(2*x$minDist+1), 1] - densStable
	dens <- dens[-(1:x$minDist), ]
	
	if(background && x$weight < 1){
		bg <- indNuc(meanDist=x$meanDist, weight=1-x$weight, length=nrow(dens))
		dens[, 1] <- dens[,1] + bg[,1]
	}
	
	apply(dens, 2, rev)
}

featureDensity.StableFeature <- function(x, stable=stableDens, background=FALSE, ...){
	dens <- stable(ceiling(-0.5*x$length):floor(0.5*x$length), shift=x$shift, ratio=x$ratio, 
			weight=x$weight, stability=x$stability)
	dens <- cbind(dens, rep(x$weight, length(dens)))
	if(background){
		dens[,1] <- dens[,1] + indNuc(meanDist=x$meanDist, weight=1-x$weight, length=nrow(dens))[,1]
	}
	dens
}

featureDensity.NFRFeature <- function(x, background=FALSE, ...){
	dens <- noNuc(length=x$length, weight=x$weight)
	
	if(background && x$weight < 1){
		bg <- indNuc(meanDist=max(x$length, 200), weight=1-x$weight, length=x$length)
		dens[, 1] <- dens[,1] + bg[,1]
	}
	dens
}

featureDensity.FuzzyFeature <- function(x, ...){
	indNuc(meanDist=x$meanDist, weight=x$weight, length=x$length)
}


## take a list of (reconciled) features and compute combined nucleosome density
## features: list of features
## length: total length of resulting density vector (may be missing)
## featureBgr: logical indicating whether feature specific background should be added
feat2dens <- function(features, length, featureBgr=TRUE, ...){
	## get densities for each feature
	featureDens <- lapply(features, featureDensity, background=featureBgr, ...)
	
	## combine everything
	if(missing(length))
		length <- sum(sapply(features, "[[", "length")) - 
				sum(sapply(features, "[[", "overlap"), na.rm=TRUE) + 
				features[[1]]$start - 1
	dens <- matrix(0, ncol=2, nrow=length)
	
	dens[features[[1]]$start:(features[[1]]$length + features[[1]]$start - 1), ] <- 
			if(is.matrix(featureDens[[1]]))featureDens[[1]] else cbind(featureDens[[1]], rep(1, length(featureDens[[1]])))
	last <- features[[1]]$length + features[[1]]$start - 1
	if(length(features) > 1) 
		for(i in 2:length(features)){
			start <- last-features[[i-1]]$overlap + 1
			dens[start:(start + features[[i]]$length -1), ] <- 
					joinRegion(dens[start:last, ], 
							if(is.matrix(featureDens[[i]]))featureDens[[i]] 
							else cbind(featureDens[[i]], rep(1, length(featureDens[[i]]))), 
					features[[i-1]]$overlap, features[[i-1]]$overlapWeight)
			last <- start + features[[i]]$length -1
		}

	dens[dens[,1] < 0 , 1] <- 0
	dens[,1]
}




######################## SECTION: generating nucleosome density ################################################

## nucleosome density for a single stable nucleosome, centred at 0
## x: position at which the density should be evaluated
## shift: size of shift for alternative nucleosome positions (to the left and right of central position)
## ratio: ratio of alternative to central position
## weight: weight of nucleosome (max 1). The density will be scaled to this value
## stability: a measure of how stable the nucleosome is (in [0, 1])
stableDens <- function(x, shift=10, ratio=1, weight=1, stability=1){
	sd <- ifelse(shift==0, 2.5, shift/4)/stability
	if(shift == 0) ratio <- 0
	min(weight, 1)*(0.5*ratio*dnorm(x, mean=-shift, sd=sd) + 
				dnorm(x, sd=sd) + 
				0.5*ratio*dnorm(x, mean=shift, sd=sd))/(ratio + 1) 
}

## distribution of distances between nucleosomes (centre to centre)
distDens <- function(x, minDist=175, varDist=337.5, meanDist=200){
	meanDist <- meanDist - minDist
	dnbinom(x-minDist, size=meanDist^2/(varDist-meanDist), mu=meanDist)
}

## distribution of fragment length
## Note that this gives the distribution of the variable part of the fragment length,
## ie it does not include the length of the binding site
## x: distance from nucleosome centre
## meanLength: mean fragment length
## minLength: minimum fragment length
## maxLength: maximum fragment length
fragDens <- function(x, minLength, maxLength, meanLength, bind){
	minLength <- minLength - bind
	maxLength <- maxLength - bind
	meanLength <- meanLength - bind
	x <- x - bind
	ifelse(x >= minLength & x <= maxLength, dnbinom(x, size=2, mu=meanLength)/(pnbinom(maxLength, size=2, mu=meanLength) - 
						pnbinom(minLength-1, size=2, mu=meanLength)), 0)
}

## Distribution of binding site within fragment
## x: position in fragment. Either absolute distance from start or 
##    fraction left of binding site centre (if fragLength is not missing)
## alpha, beta: parameters of beta distribution
## note that this assumes a symmetric distribution
bindLocDens <- function(x, fragLength){
	if(!missing(fragLength)) x <- x/fragLength 
	dbeta(x, shape1=2, shape2=2)
}

## compute nucleosome density for phased nucleosomes
## stable: function giving the density for a stable nucleosome centred at 0
## dist: function giving the distribution of distances between nucleosomes (centre to centre)
## minDist: minimum distance between nucleosomes
## TODO: handle additional arguments to stable and dist more gracefully
phaseNuc <- function(stable, dist, minDist=175, length=2000, meanDist=200, 
		varDist=(meanDist-minDist)+(meanDist-minDist)^2/2, shift=10, ratio=1, weight=1, stability=1){
	if(is.null(meanDist)) meanDist <- 200 ## for featureDensity.StableFeature
	nnuc <- floor(length/minDist)-1
	dens <- numeric(2*length)
	dens <- stable(ceiling(-0.5*length):floor(1.5*length), shift=shift, ratio=ratio, weight=weight, stability=stability)
	
	densLength <- length(dens)
	distVec <- dist(0:(densLength-1), minDist=minDist, varDist=varDist, meanDist=meanDist)
	padding <- nextn(densLength) - densLength
	distFFT <- fft(c(distVec,rep(0, padding)))
	d <- c(dens, rep(0, padding))
	for(i in 1:nnuc){
		densVec <- d
		d <- Re(fft(fft(densVec) * distFFT, inverse=TRUE)/(length(densVec)))
		dens[(i*minDist+1):densLength] <- dens[(i*minDist+1):densLength] + d[(i*minDist+1):densLength]
	}
	cbind(dens[(floor(0.5 * length)-minDist+1):(floor(1.5*length)-minDist)], rep(weight, length))
}

## nucleosome density for independent nucleosomes (ie, no phasing)
## meanDist: average distance between nucleosomes
## length: length of window
indNuc <- function(meanDist=200, length=2000, weight=1){
	cbind(weight * rep(1/meanDist, length), rep(weight, length))
}

## nucleosome free region
noNuc <- function(length, weight=1){
	cbind(rep(0, length), rep(weight, length))
}


## takes to vectors (of nucleosome densities) and stiches them together
## left, right: the two regions to be combined
## overlap: overlap between regions Defaults to 1/4 of region length
## overlapWeights: function to be used to calculate weights in overlapping region
joinRegion <- function(left, right, overlap, overlapWeights){
	if(missing(overlap)) overlap <- min(floor(length(left)/4), floor(length(right)/4))
	if(overlap == 0) return(right)
	weight <- list(left[, 2], right[, 2])
	left <- left[, 1]
	right <- right[, 1]
	
	drop <- cbind(dnorm(seq(0,5, length.out=overlap))/dnorm(0)*weight[[1]][(length(left)-overlap+1):length(left)], 
			rev(dnorm(seq(0,5, length.out=overlap))/dnorm(0)*weight[[2]][1:overlap]))
	drop <- drop / rowSums(drop)
	
	result <- numeric(length(left)+length(right)- overlap)
	result[1:(length(left)-overlap)] <- left[1:(length(left)-overlap)]
	result[(length(left)-overlap + 1):length(left)] <- drop[,1] * left[(length(left)-overlap+1):length(left)] + 
			drop[,2] * right[1:overlap]
	result[(length(left) + 1):length(result)] <- right[(overlap+1):length(right)] 

	weightNew <- numeric(length(result))
	weightNew[1:(length(left)-overlap)] <- weight[[1]][1:(length(left)-overlap)]
	weightNew[(length(left) + 1):length(result)] <- weight[[2]][(overlap+1):length(right)]
	if(!missing(overlapWeights) && !is.null(overlapWeights)) 
		weightNew[(length(left)-overlap + 1):length(left)] <- mapply(overlapWeights, 
				drop[,1] * weight[[1]][(length(left)-overlap + 1):length(left)], 
				drop[,2] * weight[[2]][1:overlap])
	else{
		weightNew[(length(left)-overlap + 1):length(left)] <- weight[[1]][(length(left)-overlap + 1):length(left)]
		if(!isTRUE(all.equal(weight[[1]][(length(left)-overlap + 1):length(left)], weight[[2]][1:overlap]))){
			warning("Unreconciled weights in overlapping region")
		}
	}
			
	rescale <- weightNew > 1
	if(any(rescale)){
		result[rescale] <- result[rescale]/weightNew[rescale]
		weightNew[rescale] <- 1
	}
	
	cbind(result, weightNew)
}


################################### SECTION: generating reads ########################################################

## convert binding site density into read density
## bindDens: vector of nucleosome densities
## fragment: function giving the fragment length distribution
## nfrag: number of fragments to simulate to generate read distribution
## ...: further arguments to fragment
bindDens2readDens <- function(bindDens, fragment, nfrag=1e5, bind=147, minLength=150, maxLength=180, ...){
	fragSample <- sample((minLength:maxLength)-bind, nfrag, replace=TRUE,
			prob=fragment((minLength:maxLength), minLength=minLength, maxLength=maxLength, bind=bind, ...))
	## get binding site distribution
	step <- 1/(maxLength-bind+1)
	locDist <- bindLocDens(seq(0,1, by=step))
	## distribution of reads on one strand (assumed to be identical on the other)
	readSample <- round(fragSample * sample(seq(0, 1, by=step), nfrag, prob=locDist, replace=TRUE))
	readDist <- hist(readSample, breaks=-1:(maxLength-bind + 1)+0.5, plot=FALSE)$density
	readDist <- readDist/sum(readDist)

	readKernel <- c(rep(0, floor(bind * 0.5)),readDist)
	readKernel <- c(readKernel, rep(0, nextn(length(readKernel)) - length(readKernel)))
	n <- length(bindDens)
	bindDens <- c(bindDens, rep(0, nextn(n + length(readKernel) - 1) - (n + length(readKernel) - 1)))
	idx <- list((length(readDist)+floor(bind * 0.5) + 1):(n + length(readDist)+floor(bind * 0.5)), 1:n)
	readDens <- cbind(convolve(bindDens, readKernel, type="open")[idx[[1]]],	
			convolve(bindDens, rev(readKernel), type="open")[idx[[2]]])
	readDens <- apply(readDens, 2, function(x) ifelse(x < 0, 0, x))
	
	## TODO: ensure reads are not generated closer to the end of the chromosome than the insert size allows
	
	readDens
}

## sample read start positions from read density
sampleReads <- function(readDens, nreads=6e6, strandProb=c(0.5, 0.5)){
	names(strandProb) <- c("fwd", "rev")
	colnames(readDens) <- c("fwd", "rev")
	readPos <- list(fwd=numeric(), rev=numeric())
		
	strand <- table(sample(names(strandProb), nreads, prob=strandProb, replace=TRUE))
	if(is.na(strand["fwd"])) strand["fwd"] <- 0
	if(is.na(strand["rev"])) strand["rev"] <- 0
	readPos$fwd <- sample(nrow(readDens), strand["fwd"], prob=readDens[, "fwd"], replace=TRUE)
	readPos$rev <- sample(nrow(readDens), strand["rev"], prob=readDens[, "rev"], replace=TRUE)
	
	readPos
}

## extract read qualities from an ShortReadQ object or from file
## reads: object or file name
## minLength: minimum length of reads for which quality is extracted
## returns list with read quality scores
extractQuality <- function(reads, minLength=25, dir, type=c("Illumina", "Sanger", "Solexa")){
	if(!is(reads, "ShortReadQ") && !missing(dir)) reads <- ShortRead::readFastq(dir, reads)
	type <- match.arg(type[1], c("Illumina", "Sanger", "Solexa"))
	offset <- c(64, 33, 59)
	names(offset) <- c("Illumina", "Sanger", "Solexa")
		
	quality <- quality(reads)
	length <- numeric(length(quality))
	for(i in 1:length(quality)){
		length[i] <- length(quality[[i]])
	}
	quality <- quality[length >= minLength]
	
	## convert quality scores to probabilities
	numericQuality <- vector(length(quality), mode="list")
	for(i in 1:length(quality)){
		numericQuality[[i]] <- as.numeric(sapply(as.character(quality[[i]]), charToRaw)) - offset[type]
	}
	## Phred score
	if(type %in% c("Illumina", "Sanger"))
		numericQuality <- lapply(numericQuality, function(x) 10^(x/-10))
	## Solexa score
	if(type %in% c("Solexa")){
		numericQuality <- lapply(numericQuality, function(x) 10^(x/-10))
		numericQuality <- lapply(numericQuality, function(x) x/(x+1))
	}
	
	numericQuality	
}

decodeQuality <- function(quality, type=c("Illumina", "Sanger", "Solexa")){
	type <- match.arg(type[1], c("Illumina", "Sanger", "Solexa"))
	offset <- c(64, 33, 59)
	names(offset) <- c("Illumina", "Sanger", "Solexa")
	
	numericQuality <- as.numeric(charToRaw(as.character(quality))) - offset[type]
	
	## Phred score
	if(type %in% c("Illumina", "Sanger"))
		numericQuality <-  10^(numericQuality/-10)
	## Solexa score
	if(type %in% c("Solexa")){
		numericQuality <- 10^(numericQuality/-10)
		numericQuality <- numericQuality/(numericQuality+1)
	}
	
	numericQuality
}

## convert error probability into integer code
encodeQuality <- function(quality, type=c("Illumina", "Sanger", "Solexa")){
	type <- match.arg(type[1], c("Illumina", "Sanger", "Solexa"))
	offset <- c(64, 33, 59)
	names(offset) <- c("Illumina", "Sanger", "Solexa")
	
	## Phred score
	if(type %in% c("Illumina", "Sanger"))
		quality <- -10 * log(quality, 10)
	## Solexa score
	if(type %in% "Solexa")
		quality <- -10 * log(quality/(1-quality), 10)
	
	rawToChar(as.raw(quality + offset[type]))
}

## determine quality scores for reads
## read: read sequence
## qualities: list of observed read qualities
## returns a read quality vector of same length as read
readQualitySample <- function(read, qualities, checkLength=TRUE, ...){
	if(checkLength) idx <- sample(which(IRanges::width(qualities) >= length(read)), 1)
	else idx <- sample(length(qualities), 1)
	XVector::subseq(qualities[[idx]], 1, length(read))
}

## probability of sequencing error producing a certain nucleotide given the correct nucleotide 
defaultErrorProb <- function(){
	prob <- list(A=c(0, 0.14, 0.05, 0.05), C=c(0.13, 0, 0.02, 0.04), G=c(0.04, 0.08, 0, 0.12), T=c(0.08, 0.15, 0.09, 0))
	prob <- lapply(prob, "names<-", c("A", "C", "G", "T"))
	prob
}

## introduce errors into read sequence based on qualities
## read: sequence of read
## qual: error probability
## prob: list with substitution probabilities for each symbol in alphabet
readError <- function(read, qual, alphabet=c("A", "C", "G", "T"), prob=defaultErrorProb(), ...){
	## replace everything not in alphabet wit first letter in alphabet
	## (with default settings this converts 'NNNNN' into 'AAAAA')
	## TODO: reduce quality for replaced nucleotides
	read <- gsub(paste("[^", paste(alphabet, collapse="", sep=""),"]", sep=""), alphabet[1], read)
	## determine location of error
	errorPos <- runif(length(qual)) < qual
	if(any(errorPos)){
		for(i in which(errorPos)){
			transProb <- prob[[substr(read, i, i)]]
			substr(read, i, i) <- sample(names(transProb), 1, prob=transProb)
		}
	}
	read
}

## get read sequence for single read from reference
readSequence <- function(readPos, sequence, strand, readLen=36){
	## get true read sequence
	if(strand == -1) readPos <- readPos - readLen + 1
	seq <- XVector::subseq(sequence, readPos, width=readLen)
	if(strand == -1) seq <- Biostrings::reverseComplement(seq)
	
	seq	
}

## write files in Illumins FASTQ format
writeIllumina <- function(reads, quality, file="", nameBase="HWSIM-EAS000R", append=FALSE){
	nlanes <- min(ceiling(sum(sapply(reads, length))/4e6), 8)
	ntiles <- 100
	clusterRange <- c(49, 2000)
	
	maxReads <- 0.5 * nlanes * ntiles * diff(clusterRange) 
	if( maxReads < sum(sapply(reads, length))){
		warning("Only first ", maxReads, " reads are used.")
		if(length(reads[[1]]) > maxReads/2) reads[[1]] <- reads[[1]][1:floor(maxReads/2)]
		reads[[2]] <- reads[[2]][1:min(maxReads-length(reads[[1]]), length(reads[[2]]))]
	}
	
	names <- lapply(1:2, function(i) vector(mode="character", length(reads[[i]])))
	for(i in 1:2){
		for(j in 1:length(reads[[i]])){
			newName <- ""
			while(newName == "" || newName %in% names[[1]] || newName %in% names[[2]]){
				newName <- paste(nameBase, sample(1:nlanes, 1), sample(1:ntiles, 1), 
						sample(clusterRange[1]:clusterRange[2],1 ), sample(clusterRange[1]:clusterRange[2],1 ), 
						sep=":")
			}
			names[[i]][j] <- newName
			cat("@", names[[i]][j], "\n", sep="", file=file, append=append || i > 1 || j > 1)
			cat(as.character(reads[[i]][j]), "\n", sep="", file=file, append=TRUE)
			cat("+", names[[i]][[j]], "\n", sep="", file=file, append=TRUE)
			cat(encodeQuality(quality[[i]][j], type="Illumina"), "\n", sep="", file=file, append=TRUE)
		}
	}
	invisible(names)
}

writeFASTQ <- function(read, quality, name, file, append=FALSE){
	if(file != "" && !append) file.create(file)
	mapply(function(r, q, n){
				cat("@", n, "\n", sep="", file=file, append=TRUE)
				cat(as.character(r), "\n", sep="", file=file, append=TRUE)
				cat("+", n, "\n", sep="", file=file, append=TRUE)
				cat(as.character(q), "\n", sep="", file=file, append=TRUE)
			}, read, quality, name)
	invisible(NULL)
}

## Convert read positions for a single chromosome (both strands) into read sequences + qualities and
## write them to file
pos2fastq <- function(readPos, names, quality, sequence, qualityFun, errorFun, readLen=36, file, 
		qualityType=c("Illumina", "Sanger", "Solexa"), ...){

	if(file != "" && !is(file, "connection") && !file.exists(file)) file <- file(file, open="w", blocking=FALSE)
	replaceProb <- if(is.null(list(...)$prob)) defaultErrorProb() else match.fun(list(...)$prob)
	qualityFun <- match.fun(qualityFun)
	errorFun <- match.fun(errorFun)
	for(i in 1:2){
		for(j in 1:length(readPos[[i]])){
			## get (true) sequence
			readSeq <- readSequence(readPos[[i]][j], sequence, strand=ifelse(i==1, 1, -1), readLen=readLen)
			## get quality
			readQual <- qualityFun(readSeq, quality, ...)
			## introduce sequencing errors
			readSeq <- errorFun(readSeq, decodeQuality(readQual, type=qualityType), prob=replaceProb)
			## write to file
			writeFASTQ(as.character(readSeq), as.character(readQual), names[[i]][j], file=file, append=TRUE)
		}
	}
		
	invisible(sum(sapply(readPos, length)))
}

writeReads <- function(readPos, readNames, sequence, quality, file, ...){
	if(is(quality, "connection") || (is.character(quality) && file.exists(quality))) 
		quality <- ShortRead::readFastq(quality)
	if(is(quality, "ShortReadQ")) quality <- quality(quality(quality))
	
	fileName <- file
	file <- file(file, open="w", blocking=FALSE)
	for(i in 1:length(sequence))
		pos2fastq(readPos=readPos[[i]], names=readNames[[i]], sequence[[i]], quality=quality, file=file, ...)
	
	close(file)
	fileName
}

## create random (Solexa style) read names
solexaNames <- function(n, nameBase="HWSIM-EAS000R"){
	nlanes <- min(ceiling(n/4e6), 8)
	ntiles <- 100
	clusterRange <- c(49, 2000)
	
	maxReads <- 0.5 * nlanes * ntiles * diff(clusterRange)^2
	if( maxReads < n){
		warning("Only ", maxReads, " read names will be generated.")
		n <- maxReads
	}
	
	names <- character(n)
	clusters <- clusterRange[1]:clusterRange[2]
	lanes <- 1:nlanes
	tiles <- 1:ntiles
	for(i in 1:n){
		newName <- ""
		while(newName == "" || newName %in% names){
			newName <- paste(nameBase, sample(lanes, 1), sample(tiles, 1), 
					sample(clusters, 1), sample(clusters, 1), 
					sep=":")
		}
		names[i] <- newName
	}
	
	names
}

## create simple read names
simpleNames <- function(n, nameBase="read"){
	paste(nameBase, 1:n, sep="_")	
}

################################# SECTION: Main driver ########################################
## This function acts as driver for the simulation. It takes all required arguments and passes
## them on to the functions for the various stages of the simulation.
## The simulation consists of a number of stages:
## 1. generate feature sequence: (genome) sequence length -> feature sequence (list)
## 2. compute nucleosome density: feature sequence [, genome length] -> nucleosome density (vector) 
## 3. compute read density: nucleosome density -> read density (matrix)
## 4. sample read start sites: read density -> read positions (list)
## 5. create read names: number of reads -> unique names
## 6. obtain read sequence and quality: read positions, genome sequence, [qualities] -> output file 

## nreads: number of reads to generate
## genome: an object of class 'DNAStringSet' or the name of a fasta file containing the genome sequence
## file: base of output file names
## functions: named list of functions to use for various stages, expected names are:
##            'features', 'nucDens', 'readDens', 'sampleReads', 'readSequence', 'formatWriter'
## control: named list of arguments to be passed to simulation functions (one list per function)
## verbose: logical indicating whether progress messages should be printed
## load: logical indicating whether an attempt should be made to load previous results
simChIP <- function(nreads, genome, file, functions=defaultFunctions(), 
		control=defaultControl(), verbose=TRUE, load=FALSE){
	## collect timing information
	timeStart <- Sys.time()
	
	## check that all functions are present
	funName <- c('features', 'bindDensity', 'readDensity', 'sampleReads', 'readNames', 'readSequence')
	if(any(!(names(functions) %in% funName))) 
		warning("Functions ", names(functions)[!(names(functions) %in% funName)]," will not be used.")
	if(!all(funName %in% names(functions))) stop("Functions for ",funName[!funName %in% names(functions)],
				" are missing.")
	
	## ensure the genome sequence is available
	if(missing(genome)) stop("Argument \"genome\" is missing, with no default. Please provide a reference sequence.")
	if(is(genome, "DNAString")) genome <-  DNAStringSet(genome)
	if(!is(genome, "DNAStringSet")) 
		genome <-  readDNAStringSet(genome, format="fasta")
	
	## ensure we know how many reads to generate
	if(missing(nreads)) stop("Argument \"nreads\" is missing, with no default. Please provide the number of ",
				"sequence reads to simulate.")
	
	## if no file name is provided we just write to the console (and don't save intermediate results)
	if(missing(file)) file <- ""

	## ensure we have legal function names or objects
	functions <- lapply(functions, function(f) if(is.function(f) || is.name(f)) f else as.name(f))
	
	## check whether there are results from a previous run that we can re-use
	stageComplete <- logical(3)
	if(load){
		if(file.access(paste(file, "features.rdata", sep="_"), mode=4) == 0) stageComplete[1] <- TRUE
		if(file.access(paste(file, "bindDensity.rdata", sep="_"), mode=4) == 0) stageComplete[2] <- TRUE
		if(file.access(paste(file, "readDensity.rdata", sep="_"), mode=4) == 0) stageComplete[3] <- TRUE
	}
	
	## check for existing result files
	resultExt <- 1
	if(nchar(file) > 0){
		resultExt <- length(Sys.glob(paste(file, "*", "fastq.txt", sep="_"))) + 1
	}
	
	## call function to generate feature sequence. We create one sequence per chromosome
	timeStage <- Sys.time()
	features <- NULL
	if(load && stageComplete[1] && !any(stageComplete[2:3])){
		warn <- getOption("warn")
		options(warn=-1)
		loaded <- try(load(paste(file, "features.rdata", sep="_")), silent=TRUE)
		if(is(loaded, "try-error"))
			stageComplete[1] <- FALSE
		if(stageComplete[1] && verbose) cat("Features loaded.")
		options(warn=warn)
	}
	if(!stageComplete[1]){
		if(verbose) cat("Generating features...")
		features <- lapply(IRanges::width(genome), function(w){
					control[[funName[1]]]$length <- w
					do.call(functions[[funName[1]]], control[[funName[1]]])
				} )
		if(nchar(file) > 0) save(features, file=paste(file, "features.rdata", sep="_"))
		stageComplete[1] <- TRUE
	}
	timeEnd <- Sys.time()
	if(verbose && !is.null(features)){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	} 
	
	## compute nucleosome density
	timeStage <- Sys.time()
	bindDensity <- NULL
	if(load && stageComplete[2] && !stageComplete[3]){
		warn <- getOption("warn")
		options(warn=-1)
		loaded <- try(load(paste(file, "bindDensity.rdata", sep="_")), silent=TRUE)
		if(is(loaded, "try-error")) 
			stageComplete[2] <- FALSE
		if(stageComplete[2] && verbose) cat("Binding site density loaded. ")
		options(warn=warn)
	}
	if(!stageComplete[2]){
		if(verbose) cat("Computing binding site density... ")
		bindDensity <- mapply(function(f, w){
					control[[funName[2]]]$features <- f
					control[[funName[2]]]$length <- w
					eval(as.call(c(functions[[funName[2]]], control[[funName[2]]])))
				}, features, IRanges::width(genome), SIMPLIFY=FALSE)
		if(nchar(file) > 0) save(bindDensity, file=paste(file, "bindDensity.rdata", sep="_"))
		stageComplete[2] <- TRUE
	}
	if(nchar(file) > 0) features <- paste(file, "features.rdata", sep="_")
	timeEnd <- Sys.time()
	if(verbose && !is.null(bindDensity)){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	}
	
	## compute read density
	timeStage <- Sys.time()
	readDensity <- NULL
	if(load && stageComplete[3]){
		warn <- getOption("warn")
		options(warn=-1)
		loaded <- try(load(paste(file, "readDensity.rdata", sep="_")), silent=TRUE)
		if(is(loaded, "try-error")){
			stageComplete[3] <- FALSE
		}
		if(stageComplete[3] && verbose) cat("Read density loaded. ")
		options(warn=warn)
	}
	if(!stageComplete[3]){
		if(verbose) cat("Computing read density... ")
		readDensity <- lapply(bindDensity, 
				function(d){
					control[[funName[3]]]$bindDens <- d
					eval(as.call(c(functions[[funName[3]]], control[[funName[3]]])))
				})
		if(nchar(file) > 0) save(readDensity, file=paste(file, "readDensity.rdata", sep="_"))
		stageComplete[3] <- TRUE
	}
	if(nchar(file) > 0) bindDensity <- paste(file, "bindDensity.rdata", sep="_")
	timeEnd <- Sys.time()
	if(verbose && !is.null(readDensity)){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	}
	
	## sample reads
	timeStage <- Sys.time()
	if(verbose) cat("Sampling reads... ")
	readPosition <- mapply( 
			function(d, w){
				control[[funName[4]]]$readDens <- d
				control[[funName[4]]]$nreads <- round(nreads * w / sum(IRanges::width(genome)))
				eval(as.call(c(functions[[funName[4]]], control[[funName[4]]])))
			}, readDensity, IRanges::width(genome), SIMPLIFY=FALSE)
	nreads <- sum(sapply(readPosition, sapply, length))
	if(nchar(file) > 0){ 
		save(readPosition, file=paste(file, resultExt, "readPosition.rdata", sep="_"))
		readDensity <- paste(file, "readDensity.rdata", sep="_")
	}
	timeEnd <- Sys.time()
	if(verbose){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	}
	
	## generate read names
	timeStage <- Sys.time()
	if(verbose) cat("Generating read names... ")
	control[[funName[5]]]$n <- nreads
	readNames <- eval(as.call(c(functions[[funName[5]]], control[[funName[5]]])))
	## reshape to match readPosition
	readNames <- relist(readNames, readPosition)
	if(nchar(file) > 0) save(readNames, file=paste(file, resultExt, "readNames.rdata", sep="_"))
	timeEnd <- Sys.time()
	if(verbose){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	}
		
	## get read sequence and quality
	timeStage <- Sys.time()
	if(verbose) cat("Determining read sequence and quality... ")
	control[[funName[6]]]$readPos <- readPosition
	control[[funName[6]]]$readNames <- readNames
	control[[funName[6]]]$file <- paste(file, resultExt, "fastq.txt", sep="_")
	control[[funName[6]]]$sequence <- genome
	readSequence <- do.call(functions[[funName[6]]], control[[funName[6]]])
	timeEnd <- Sys.time()
	if(verbose){
		thisStage <- difftime(timeEnd, timeStage)
		cat("(", thisStage, " ", units(thisStage), ")\n", sep="")
	}
	if(nchar(file) > 0) readNames <- paste(file, resultExt, "readNames.rdata", sep="_")

	timeEnd <- Sys.time()
	if(verbose){
		total <- difftime(timeEnd, timeStart)
		cat("Total time:", total, units(total), "\n")
	}
	
	if(nchar(file) > 0) readPosition <- paste(file, resultExt, "readPosition.rdata", sep="_")
	list(features=features, bindDensity=bindDensity, readDensity=readDensity, readPosition=readPosition, 
			readSequence=readSequence, readNames=readNames)	
}


defaultFunctions <- function(){
	list(features=placeFeatures, bindDensity=feat2dens, readDensity=bindDens2readDens, sampleReads=sampleReads, 
			readSequence=writeReads, readNames=simpleNames)
}

defaultControl <- function(features=list(), bindDensity=list(), readDensity=list(), 
		readNames=list(), readSequence=list()){
	fragment <- readDensity$fragment
	meanLength <- readDensity$meanLength
	readDensity <- readDensity[!(names(readDensity) %in% c("fragment", "meanLength"))]
	readDensity <- c(list(fragment=if(is.null(fragment)) fragDens else fragment, 
					meanLength=if(is.null(meanLength)) 160 else meanLength), readDensity)
	
	qualityFun <- readSequence$qualityFun
	errorFun <- readSequence$errorFun
	readLen <- readSequence$readLen
	readSequence <- readSequence[!(names(readSequence) %in% c("qualityFun", "errorFun", "readLen"))]
	readSequence <- c(list(qualityFun=if(is.null(qualityFun))readQualitySample else qualityFun, 
					errorFun=if(is.null(errorFun))readError else errorFun, 
					readLen=if(is.null(readLen)) 36 else readLen), readSequence)
	
	list(features=features, bindDensity=bindDensity, readDensity=readDensity, 
			readNames=readNames, readSequence=readSequence)
}

## randomly choose records from a file
## file: input file
## n: number of records to select
## nrec: total number of records in file
## recLen: length of record (in lines)
## output: output file name
sampleFromFile <- function(file, n, nrec, recLen=4, skip=0, output){
	if(missing(nrec)) nrec <- (ShortRead::countLines(file)-skip)/recLen
	file <- file(file)
	open(file)
	readLines(file, skip)
	file.create(output)
	output <- file(output) 
	open(output, open="w")
	
	## select records
	selected <- sort(sample(nrec, n, replace=FALSE), decreasing=FALSE)
	idx <- 1
	for(i in 1:selected[length(selected)]){
		record <- readLines(file, recLen)
		if(i == selected[idx]){
			writeLines(record, output)
			idx <- idx + 1
		}
	}
	
	close(file)
	close(output)
}

#################### SECTION: print methods ############################
print.SimulatedExperiment <- function(x, ...){
	classes <- sort(unique(sapply(x, function(y) class(y)[1])))
	cat("Object of class ", class(x)[1], if(inherits(x, "ReconciledSimulatedExperiment")) " (reconciled)", 
			" with ", length(x), " feature", sep="")
	if(length(classes) > 1) cat("s")
	cat(".\nFeature class")
	if(length(classes) > 1) cat("es")
	cat(": ",classes,"\n")
	invisible(x)
}

print.SimulatedFeature <- function(x, ...){
	cat("Object of class", class(x)[1], if(inherits(x, "ReconciledFeature")) "(reconciled)", "\n")
	for(i in 1:length(x)){
		cat(names(x)[i], ": ", sep="")
		if(!is(x[[i]], "function")) cat(x[[i]], "\n")
		else str(x[[i]])
	}
	invisible(x)
}


#################### SECTION: summary methods ############################
summary.SimulatedExperiment <- function(object, ...){
	classes <- sapply(object, function(x) class(x)[1])
	cat("Object of class ", class(object)[1], if(inherits(object, "ReconciledSimulatedExperiment")) " (reconciled)", 
			" with ", length(object), " feature", sep="")
	if(length(classes) > 1) cat("s")
	cat("\n")
	x <-table(classes)
	for(i in 1:length(x)) cat(names(x)[i], ": ", x[[i]], "\n", sep="")
	invisible(x)
}
