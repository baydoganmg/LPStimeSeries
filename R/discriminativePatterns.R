discriminativePatterns <- function(object, x, classes, n=5, plot=TRUE,
                                   orient=c(2,2), palette=NULL) {

  if (!inherits(object, "learnPattern"))
    stop("object not of class learnPattern")
  if (is.null(object$forest))
    stop("No forest component in the object")

  if (!is.matrix(x)) {
    if (length(x) > 0) {
      x <- t(as.matrix(x))
    } else {
      stop("data (x) has 0 rows")
    }
  }

  classes <- as.factor(classes)
  if (length(classes) != nrow(x))
    stop("length of classes must equal the number of time series (rows of x)")
  if (nlevels(classes) < 2)
    warning("fewer than 2 classes; all discriminability scores will be zero")

  n <- max(1L, as.integer(n))

  ## Compute representation (terminal node frequencies)
  rep <- predict(object, x, nodes=TRUE)

  ## Per-class mean frequency for each terminal node
  class_levels <- levels(classes)
  nclasses <- nlevels(classes)
  class_sizes <- as.numeric(table(classes))
  class_means <- t(sapply(class_levels, function(cl)
    colMeans(rep[classes == cl, , drop=FALSE])))
  overall_means <- colMeans(rep)

  ## Discriminability score: weighted between-class variance / mean
  deviations <- sweep(class_means, 2, overall_means)^2
  weighted_devs <- sweep(deviations, 1, class_sizes, "*")
  scores <- colSums(weighted_devs) / pmax(overall_means, .Machine$double.eps)
  scores[overall_means == 0] <- 0

  ## Rank and select top-n
  ranking <- order(scores, decreasing=TRUE)
  top_n <- ranking[seq_len(min(n, length(ranking)))]

  ## Plot top patterns colored by class
  if (plot && nlevels(classes) >= 2) {
    if (is.null(palette)) {
      palette <- if (nclasses < 12)
        brewer.pal(max(3, nclasses), "Set1")[1:nclasses] else rainbow(nclasses)
    }

    nofterminals <- c(1, apply(object$forest$nodestatus, 2,
                               function(ns) sum(ns == -1)))
    nofterminals <- cumsum(nofterminals)
    mdim <- ncol(x)
    ntest <- nrow(x)
    xt <- t(data.matrix(x))

    nplots <- min(length(top_n), prod(orient))
    op <- par(mfrow=orient, mar=c(4, 3, 3, 1))
    on.exit(par(op))

    ## Faint colors for background series and shading per class
    palette_bg <- adjustcolor(palette, alpha.f=0.15)
    palette_shade <- adjustcolor(palette, alpha.f=0.08)

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

    for (p in seq_len(nplots)) {
      tid <- top_n[p]

      ## Map global terminal node ID to tree + within-tree terminal
      whichtree <- findInterval(tid, nofterminals)
      terminal <- tid - nofterminals[whichtree] + 1

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
                PACKAGE = "LPStimeSeries")[c("predictpatterns", "targetpatterns")]

      ans$targetpatterns[ans$targetpatterns == -999] <- NA
      ans$predictpatterns[ans$predictpatterns == -999] <- NA
      resT <- t(matrix(ans$targetpatterns, nrow=mdim))
      resP <- t(matrix(ans$predictpatterns, nrow=mdim))

      ## Find predictor and target time windows across all series
      tgt_all <- which(colSums(!is.na(resT)) > 0)
      pred_all <- which(colSums(!is.na(resP)) > 0)

      ## y-limits from the data
      ylims <- range(x, na.rm=TRUE)
      ylims <- ylims + c(-1, 1) * 0.08 * diff(ylims)

      ## Empty plot frame
      plot(NULL, xlim=c(1, mdim), ylim=ylims,
           main=paste0("Pattern #", p, " (node ", tid,
                       ", score=", round(scores[tid], 1), ")"),
           xlab="Time", ylab="")

      ## Draw faint background series colored by class
      for (cl in seq_along(class_levels)) {
        idx <- which(classes == class_levels[cl])
        for (i in idx)
          lines(1:mdim, x[i, ], col=palette_bg[cl], lwd=0.5)
      }

      ## Shade predictor and target windows
      if (length(pred_all) > 0)
        rect(min(pred_all) - 0.5, ylims[1], max(pred_all) + 0.5, ylims[2],
             col=adjustcolor("#377EB8", alpha.f=0.08), border=NA)
      if (length(tgt_all) > 0)
        rect(min(tgt_all) - 0.5, ylims[1], max(tgt_all) + 0.5, ylims[2],
             col=adjustcolor("#E41A1C", alpha.f=0.08), border=NA)

      ## Overlay pattern segments as lines, colored by class
      for (cl in seq_along(class_levels)) {
        idx <- which(classes == class_levels[cl])
        for (i in idx) {
          tgt_pos <- which(!is.na(resT[i, ]))
          pred_pos <- which(!is.na(resP[i, ]))
          if (length(tgt_pos) > 1)
            lines_consecutive(tgt_pos, resT[i, tgt_pos], col=adjustcolor(palette[cl], alpha.f=0.5), lwd=1.5)
          if (length(pred_pos) > 1)
            lines_consecutive(pred_pos, resP[i, pred_pos], col=adjustcolor(palette[cl], alpha.f=0.5), lwd=1.5, lty=2)
        }
      }

      ## Per-class mean pattern as thick lines
      for (cl in seq_along(class_levels)) {
        idx <- which(classes == class_levels[cl])
        if (length(idx) > 1) {
          tgt_mean <- colMeans(resT[idx, , drop=FALSE], na.rm=TRUE)
          pred_mean <- colMeans(resP[idx, , drop=FALSE], na.rm=TRUE)
          tgt_ok <- which(is.finite(tgt_mean))
          pred_ok <- which(is.finite(pred_mean))
          if (length(tgt_ok) > 1)
            lines_consecutive(tgt_ok, tgt_mean[tgt_ok], col=palette[cl], lwd=3)
          if (length(pred_ok) > 1)
            lines_consecutive(pred_ok, pred_mean[pred_ok], col=palette[cl], lwd=3, lty=2)
        }
      }

      legend("topleft",
             legend=c(paste(class_levels, "(target)"),
                      paste(class_levels, "(predictor)")),
             col=rep(palette[1:nclasses], 2),
             lwd=2, lty=rep(c(1, 2), each=nclasses),
             cex=0.7, bg="white")
    }
  }

  result <- list(scores=scores, ranking=ranking, top=top_n,
                 class.means=class_means, overall.means=overall_means,
                 classes=class_levels)
  invisible(result)
}
