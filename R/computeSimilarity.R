computeSimilarity <- function(object=NULL,testseries=NULL,refseries=NULL, maxdepth=NULL, which.tree=NULL,
							sim.type=0, terminal=TRUE,testrepresentation,refrepresentation,
							nthreads=1, normalize=FALSE) {

  nthreads <- as.integer(nthreads)
  if (is.na(nthreads) || nthreads < 1) nthreads <- 1L

  if(!is.null(object)){
	  if (!inherits(object, "learnPattern"))
			stop("object not of class learnPattern")
	  if (is.null(object$forest)) stop("No forest component in the object")
	  if (is.null(refseries) && is.null(testseries)) stop("You need to provide two time series data sets for similarity computation")
	  terminal <- FALSE
  }
  
  if(terminal){ #if computation will be performed over representations
	  ntest <- nrow(testrepresentation)
	  ntrain <- nrow(refrepresentation)
	  nofterminal <- ncol(refrepresentation)
	  ans <- .C("compute_similarity", 
				as.integer(as.matrix(testrepresentation)), 
				as.integer(ntest), 
				as.integer(as.matrix(refrepresentation)),  
				as.integer(ntrain), 
				as.integer(nofterminal),
				as.integer(sim.type), 
				result=integer(ntest*ntrain))
	  return(matrix(ans$result,ntest,ntrain))
	  
  } else {	 #if computation will be performed over raw time series
  
      if(!is.matrix(testseries)){
		if(length(testseries)>0){ #single time series
			testseries <- t(as.matrix(testseries))
		}
		else{
			stop("data (testseries) has 0 rows")
		}   
    } 
      if(!is.matrix(refseries)){
		if(length(refseries)>0){ #single time series
			refseries<- t(as.matrix(refseries))
		}
		else{
			stop("data (refseries) has 0 rows")
		}   
    }     
    if(is.null(maxdepth)) maxdepth <- object$maxdepth
    if(maxdepth>object$maxdepth) {
		maxdepth <- object$maxdepth
		warning("invalid depth: reset to the maximum depth provided during training!")
    }
    if(!is.null(which.tree)){
		if(length(which.tree)==0) stop("No trees are selected!")
		usedtrees=array(0,object$ntree)
		usedtrees[which.tree]=1
	} else {
		usedtrees=array(1,object$ntree)
	}
	mdim <- ncol(refseries)
	ntree <- object$ntree
	ntrain <- nrow(refseries)
	ntest <- nrow(testseries)

	## Detect variable-length series (trailing NAs)
	if (any(is.na(refseries))) {
	  serieslens_ref <- apply(refseries, 1, function(row) {
	    na_pos <- which(is.na(row))
	    if (length(na_pos) == 0) return(length(row))
	    first_na <- min(na_pos)
	    if (any(!is.na(row[first_na:length(row)])))
	      stop("NAs in refseries must be trailing (variable-length format)")
	    return(first_na - 1L)
	  })
	  refseries[is.na(refseries)] <- 0
	} else {
	  serieslens_ref <- rep(mdim, ntrain)
	}
	if (any(is.na(testseries))) {
	  serieslens_test <- apply(testseries, 1, function(row) {
	    na_pos <- which(is.na(row))
	    if (length(na_pos) == 0) return(length(row))
	    first_na <- min(na_pos)
	    if (any(!is.na(row[first_na:length(row)])))
	      stop("NAs in testseries must be trailing (variable-length format)")
	    return(first_na - 1L)
	  })
	  testseries[is.na(testseries)] <- 0
	} else {
	  serieslens_test <- rep(mdim, ntest)
	}

	keepIndex <- c("similarity")
	x <- t(data.matrix(refseries))
	xtst <- t(data.matrix(testseries))
	ans <- .C("regForest_similarity",
			as.double(x),
			as.double(xtst),
		    as.integer(ntrain),
			as.integer(ntest),
			as.double(object$segment.length),
		    as.integer(mdim),
			as.integer(object$ntree),
			as.integer(usedtrees),
			object$forest$leftDaughter,
			object$forest$rightDaughter,
			object$forest$nodestatus,
			object$forest$nodedepth,
			object$forest$nrnodes,
			object$forest$xbestsplit,
			object$forest$bestvar,
			object$forest$splitType,
			object$forest$ndbigtree,
			as.integer(maxdepth),
			as.integer(sim.type),
			similarity = double(ntest*ntrain),
			as.integer(nthreads),
			as.integer(serieslens_ref),
			as.integer(serieslens_test),
			as.integer(normalize),
			PACKAGE = "LPStimeSeries")[keepIndex]
	}
	return(matrix(ans$similarity,ntest,ntrain))
}
