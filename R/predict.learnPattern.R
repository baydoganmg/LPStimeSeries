"predict.learnPattern" <-
    function (object, newdata, which.tree=NULL,
				nodes=TRUE, maxdepth=NULL, normalize=FALSE, ...)
{
    if (!inherits(object, "learnPattern"))
        stop("object not of class learnPattern")
    if (is.null(object$forest)) stop("No forest component in the object")
	if (!nodes&&sum(object$target.type==2)>0)  stop("Difference series are used as target segments")
	
    x <- newdata

    if (!is.numeric(x)) stop("newdata is not numeric")

    if(!is.matrix(x)){
		if(length(x)>0){ #single time series
			x <- t(as.matrix(x))
		}
		else{
			stop("data (x) has 0 rows")
		}
    }

    ## Detect variable-length series (trailing NAs)
    if (any(is.na(x))) {
      serieslens <- apply(x, 1, function(row) {
        na_pos <- which(is.na(row))
        if (length(na_pos) == 0) return(length(row))
        first_na <- min(na_pos)
        if (any(!is.na(row[first_na:length(row)])))
          stop("NAs in newdata must be trailing (variable-length format)")
        return(first_na - 1L)
      })
      if (any(serieslens == 0)) stop("Series with length 0 not permitted")
      x[is.na(x)] <- 0
    } else {
      serieslens <- rep(ncol(x), nrow(x))
    }

    if(is.null(maxdepth)) maxdepth <- object$maxdepth
    
    if(maxdepth>object$maxdepth) {
		maxdepth <- object$maxdepth
		warning("invalid depth: reset to the maximum depth provided during training!")
    }
		
    keep <- 1:nrow(x)
    rn <- rownames(x)
    if (is.null(rn)) rn <- keep

    mdim <- ncol(x)
    ntest <- nrow(x)
    
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))
    x <- t(data.matrix(x))

	if(!is.null(which.tree)){
		if(length(which.tree)==0) stop("No trees are selected!")
		usedtrees=array(0,object$ntree)
		usedtrees[which.tree]=1
	} else {
		usedtrees=array(1,object$ntree)
	}
		
	if (nodes){
		keepIndex <- c("nodeRep","lenRep")
		if(!is.null(which.tree)){
			nodexts <- double(ntest * length(which.tree) * object$forest$nrnodes)
		} else {
			nodexts <- double(ntest * object$forest$nrnodes * object$ntree )
		}

		ans <- .C("regForest_represent",
				as.double(x),
				as.integer(ntest),
				as.integer(which.tree),
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
				nodeRep = nodexts,
				lenRep = integer(1),
				as.integer(serieslens),
				as.integer(normalize),
				PACKAGE = "LPStimeSeries")[keepIndex]
				
		res=t(matrix(ans$nodeRep[1:(ans$lenRep*ntest)], nrow=ans$lenRep))
		
	} else {	

		keepIndex <- c("predicted","count")
		ans <- .C("regForest_predict",
				as.double(x),
				as.integer(ntest),
				as.integer(which.tree),
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
				as.integer(object$forest$bestvar),
				as.integer(object$forest$splitType),
				as.double(object$forest$nodepred),
				as.integer(object$forest$ndbigtree),
				as.integer(object$target),
				as.integer(maxdepth),
				predicted = double(ntest * mdim),
				count = integer(ntest * mdim),
				as.integer(serieslens),
				PACKAGE = "LPStimeSeries")[keepIndex]

		ans$predicted[ans$predicted==-999]=NA
		res=list(predictions=t(matrix(ans$predicted, nrow=mdim)),
		         target.count=t(matrix(ans$count, nrow=mdim)))
	}
    res
} 
    
				   		                    
