visualizePattern <- function(object, x, which.terminal, orient=c(2,2)) {
  if (is.null(object$forest)) {
    stop("No forest component in ", deparse(substitute(object)))
  }

  if(!is.matrix(x)){
	if(length(x)>0){ #single time series
		x <- t(as.matrix(x))
	}
	else{
		stop("data (x) has 0 rows")
	}
  }

  if(is.null(which.terminal)){
	  stop("Terminal node info is not provided")
  }

  nofterminals <- c(1,apply(object$forest$nodestatus,2, function(x) sum(x==-1)))
  nofterminals <- cumsum(nofterminals)

  if(which.terminal > nofterminals[object$ntree+1]){
	  stop(paste("Total number of terminal nodes is",nofterminals[object$ntree+1],"which is less than",which.terminal))
  }

  whichtree <- findInterval(which.terminal,nofterminals)
  terminal <- which.terminal - nofterminals[whichtree] + 1
  mdim <- ncol(x)
  ntest <- nrow(x)
  xt <- t(data.matrix(x))

  keepIndex <- c("predictpatterns", "targetpatterns")
  ans <- .C("regForest_pattern",
			as.double(xt),
			as.integer(ntest),
			as.integer(whichtree),
			as.integer(terminal),
			as.double(object$segment.length),
			as.integer(mdim),
			as.integer(object$ntree),
			object$forest$leftDaughter,
			object$forest$rightDaughter,
			object$forest$nodestatus,
			object$forest$nodedepth,
			object$forest$nrnodes,
			object$forest$xbestsplit,
			object$forest$bestvar,
			object$forest$splitType,
			object$forest$ndbigtree,
			as.integer(object$maxdepth),
			as.integer(object$target),
			as.integer(object$target.type),
			predictpatterns = double(ntest * mdim),
			targetpatterns = double(ntest * mdim),
		PACKAGE = "LPStimeSeries")[keepIndex]

	ans$targetpatterns[ans$targetpatterns==-999] <- NA
	resT <- t(matrix(ans$targetpatterns,nrow=mdim))

	ans$predictpatterns[ans$predictpatterns==-999] <- NA
	res <- t(matrix(ans$predictpatterns,nrow=mdim))

	bu <- apply(res,1,function(x) sum(!is.na(x)))
	ind <- order(-bu)

    nplots <- min(prod(orient), ntest)
    op <- par(mfrow=orient, mar=c(4, 3, 3, 1))
    on.exit(par(op))

    col_target <- "#E41A1C"
    col_pred   <- "#377EB8"
    col_bg     <- adjustcolor("grey50", alpha.f=0.25)
    col_tgt_shade  <- adjustcolor(col_target, alpha.f=0.12)
    col_pred_shade <- adjustcolor(col_pred,   alpha.f=0.12)

    ## Draw lines only between consecutive time positions
    lines_consecutive <- function(pos, vals, ...) {
      if (length(pos) < 2) return()
      breaks <- which(diff(pos) > 1)
      starts <- c(1, breaks + 1)
      ends <- c(breaks, length(pos))
      for (s in seq_along(starts)) {
        seg <- starts[s]:ends[s]
        if (length(seg) > 1)
          lines(pos[seg], vals[seg], ...)
      }
    }

    for(k in 1:nplots){
      i <- ind[k]
      ylims <- range(x[i, ], na.rm=TRUE)
      ylims <- ylims + c(-1, 1) * 0.08 * diff(ylims)

      ## Empty plot frame
      plot(NULL, xlim=c(1, mdim), ylim=ylims,
           main=paste0("Pattern ", k, " (series ", i, ")"),
           xlab="Time", ylab="")

      ## Full series as faint background line
      lines(1:mdim, x[i, ], col=col_bg, lwd=1)

      ## Shade predictor window
      pred_pos <- which(!is.na(res[i, ]))
      if (length(pred_pos) > 0) {
        rect(min(pred_pos) - 0.5, ylims[1], max(pred_pos) + 0.5, ylims[2],
             col=col_pred_shade, border=NA)
        lines_consecutive(pred_pos, res[i, pred_pos], col=col_pred, lwd=2)
        points(pred_pos, res[i, pred_pos], col=col_pred, pch=17, cex=0.7)
      }

      ## Shade target window
      tgt_pos <- which(!is.na(resT[i, ]))
      if (length(tgt_pos) > 0) {
        rect(min(tgt_pos) - 0.5, ylims[1], max(tgt_pos) + 0.5, ylims[2],
             col=col_tgt_shade, border=NA)
        lines_consecutive(tgt_pos, resT[i, tgt_pos], col=col_target, lwd=2)
        points(tgt_pos, resT[i, tgt_pos], col=col_target, pch=16, cex=0.7)
      }

      legend("topleft", c("Target", "Predictor"),
             col=c(col_target, col_pred), lwd=2,
             pch=c(16, 17), cex=0.7, bg="white")
	}

	out=list(predictor=res,target=resT,tree=whichtree,terminal=terminal)
	invisible(out)
}
