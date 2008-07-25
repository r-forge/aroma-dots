fitRocLite <- function(truth, data, recall=NULL, idxs=NULL, ncuts=1200, hasNAs=TRUE, isOrdered=FALSE, ...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'data':
  if (length(data) != length(truth)) {
    throw("Argument 'data' and 'truth' are of different lengths: ", 
                               length(data), " != ", length(truth));
  }

  # Argument 'idxs':
  if (!is.null(idxs)) {
    idxs <- Arguments$getIndices(idxs);
  }


  verbose && enter(verbose, "Fitting ROC");

  truth <- as.vector(truth);
  data <- as.vector(data);

  # Remove missing values?
  if (hasNAs) {
    verbose && enter(verbose, "Removing NAs");
    ok <- (is.finite(truth) & is.finite(data));
    truth <- truth[ok];
    data <- data[ok];
    rm(ok);
    verbose && exit(verbose);
  }

  if (!isOrdered) {
    verbose && enter(verbose, "Ordering data");
    o <- order(data);
    truth <- truth[o];
    data <- data[o];
    rm(o);
    verbose && exit(verbose);
  }

  if (is.null(idxs)) {
    verbose && enter(verbose, "Identifying indices for cuts");
    n <- length(data);
    idxs <- seq(from=1, to=n, length.out=ncuts+1);
    verbose && exit(verbose);
  }

  cuts <- data[idxs];
  rm(data);

  # Turn 'truth' into (0,1) variable by re-calling?
  if (!is.null(recall)) {
    verbose && enter(verbose, "Recalling truth to (0,1)");
    verbose && cat(verbose, "Calling value: ", recall);
    truth <- (truth == recall);
    truth <- as.integer(truth);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Counting TPs and FPs in blocks");
  ncounts <- length(idxs)-1;
  counts <- matrix(as.integer(0), nrow=ncounts, ncol=2);
  colnames(counts) <- c("fpRate", "tpRate");
  for (kk in 1:ncounts) {
    ii <- idxs[kk]:(idxs[kk+1]);
    counts[kk,2] <- sum(truth[ii]);
    counts[kk,1] <- length(ii);
  }
  counts[,1] <- counts[,1] - counts[,2];
  rm(truth);

  counts <- apply(counts, MARGIN=2, FUN=cumsum);
  counts[,1] <- counts[,1] / counts[ncounts,1];
  counts[,2] <- counts[,2] / counts[ncounts,2];
  verbose && exit(verbose);


  # Return structure
  res <- list(roc=counts, cuts=cuts);

  verbose && exit(verbose);

  res;
} # fitRocLite()



############################################################################
# HISTORY:
# 2007-08-24
# o Added findTpAtFpLite() which is quite fast.
# 2007-08-20
# o Added file caching to fitRoc2().
# 2007-08-19
# o Renamed argument 'call' to 'toCall' in fitRoc().
# 2007-04-15
# o Added scanTpAtFp().
# 2007-04-14
# o Added interpolation in findTpAtFp().
# o Removed gc() from fitRoc().
# o Added findTpAtFp() to locate TP rate for given FP rate.
# 2007-03-2x
# o Added fitRoc().
# o Created.
############################################################################
