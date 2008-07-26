findTpAtFpLite <- function(truth, data, fpRate, acc=1e-3, recall=NULL, hasNAs=TRUE, isOrdered=FALSE, ..., .checkArgs=TRUE, verbose=FALSE) {
  if (.checkArgs) {
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
  
    # Argument 'fpRate':
    fpRate <- Arguments$getDouble(fpRate, range=c(0,1));
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

  # Turn 'truth' into (0,1) variable by re-calling?
  if (!is.null(recall)) {
    verbose && enter(verbose, "Recalling truth to (0,1)");
    verbose && cat(verbose, "Calling value: ", recall);
    truth <- (truth == recall);
    truth <- as.integer(truth);
    verbose && exit(verbose);
  }

  nbrOfDataPoints <- length(truth);
  totalTpCount <- as.integer(sum(truth));
  totalFpCount <- nbrOfDataPoints - totalTpCount;

  fpCount <- fpRate * totalFpCount;

  C <- as.integer(1) - truth;
  verbose && str(verbose, C);
  C <- cumsum(C);
  verbose && str(verbose, C);
  verbose && str(verbose, idx);

  idx <- match(floor(fpCount), C);
  verbose && str(verbose, idx);

  fpCountEstRange <- c(C[idx], C[idx]+1);
  fpRateEstRange <- fpCountEstRange/totalFpCount;
  d <- abs(fpRateEstRange - fpRate);
  w <- rev(d) / sum(d);
  fpRateEst <- w[1]*fpRateEstRange[1]+w[2]*fpRateEstRange[2];

  tpCountEstRange <- c(idx, idx+1) - fpCountEstRange;
  tpRateEstRange <- tpCountEstRange/totalTpCount;
  tpRateEst <- w[1]*tpRateEstRange[1]+w[2]*tpRateEstRange[2];
  verbose && exit(verbose);

  list(tpRateEst=tpRateEst, tpRateEstRange=tpRateEstRange, fpRateEst=fpRateEst, fpRateEstRange=fpRateEstRange, fpRate=fpRate, orderedIdxs=as.integer(c(idx, idx+1)), w=w);
} # findTpAtFpLite()


## IDEAS: /HB 2008-07-24
## whichLargestLessOrEqual <- function(x, threshold, ...) {
## } # whichLargestLessOrEqual()


############################################################################
# HISTORY:
# 2008-07-24
# o Created.
############################################################################
