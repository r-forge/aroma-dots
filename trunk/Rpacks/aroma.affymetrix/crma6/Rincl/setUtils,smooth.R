# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Public functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getSmoothDataSet <- function(set, h, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(2, Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Creating smoothed data set");
  verbose && cat(verbose, "Amount of smoothing: ", h);

  verbose && enter(verbose, "Extracting truth and data");
  set <- addRocData(set);
  rocData <- set$rocData;
  verbose && enter(verbose, "Calculating data checksum for file cache");
  rocDataChecksum <- getChecksum(rocData);
  verbose && exit(verbose);
  C <- getTruth(rocData, raw=TRUE);
  verbose && cat(verbose, "Truth:");
  verbose && str(verbose, C);
  M <- getData(rocData, raw=TRUE);
  verbose && cat(verbose, "Data:");
  verbose && str(verbose, M);
  rm(rocData);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting genomic positions");
  # Initial data ordered along the genome
  x <- set$positions;
  if (is.null(x))
    x <- positions;
  verbose && cat(verbose, "Positions: ");
  verbose && str(verbose, x);
  verbose && exit(verbose);

  verbose && enter(verbose, "Validating (truth, data, position) data");
  n <- length(x);
  if (nrow(M) != n) {
    throw("Number of rows in 'data' does not match the number of loci: ", nrow(M), " != ", n);
  }
  if (is.matrix(C) && nrow(C) != n) {
    throw("Number of rows in 'truth' does not match the number of loci: ", nrow(C), " != ", n);
  }
  if (is.matrix(C) && ncol(C) != ncol(M)) {
    throw("Number of columns in 'truth' and 'data' does not match: ", nrow(C), " != ", ncol(M));
  }
  verbose && exit(verbose);

  key <- list(method="getSmoothDataSet", rocDataChecksum=rocDataChecksum, x=x, h=h);
  dirs <- c("crma6");
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "Found cached results!");
    verbose && exit(verbose);
    return(res);
  }

  verbose && enter(verbose, "Ordering data along the genome");
  o <- order(x);
  x <- x[o];
  M <- M[o,,drop=FALSE];
  if (is.matrix(C))
    C <- C[o,,drop=FALSE];
  rm(o);
  verbose && cat(verbose, "Genomic positions:");
  verbose && str(verbose, x);
  verbose && cat(verbose, "Log-ratios:");
  verbose && str(verbose, M);
  verbose && cat(verbose, "Truth:");
  verbose && str(verbose, C);
  verbose && exit(verbose);


  verbose && enter(verbose, "Smoothing data along the genome");
  verbose && cat(verbose, "Amount of smoothing: ", h);

  idxs <- getBlockAverageMap(n=nrow(M), h=h);
  hApprox <- attr(idxs, "hApprox");

  x <- NULL;
  if (!is.null(x)) {
    verbose && enter(verbose, "Smoothing genome locations");
    x <- blockAvg(x, idxs);
    verbose && exit(verbose);
  }
  rm(x);

  if (is.matrix(C)) {
    verbose && enter(verbose, "Smoothing truth");
    C <- blockAvg(C, idxs);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Smoothing log-ratios");
  M <- blockAvg(M, idxs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Creating RocData object");
  rocData <- RocData(truth=C, data=M);
  rm(C, M);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating data set object");
  set <- list(
    name = sprintf("%s,h=%.3f", set$name, h),
    rocData = rocData,
    h = h, 
    hApprox = hApprox
  );
  rm(rocData);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Cache results?
  if (cache) {
    saveCache(set, key=key, dirs=dirs);
  }

  set;
} # getSmoothDataSet()


addSmoothDataSets <- function(sets, hs=2:4, subset="all", ..., verbose=FALSE) {
  # Argument 'hs':
  hs <- Arguments$getDoubles(hs, range=c(1, Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Add raw sets for subset, if missing
  sets <- addSubsetDataSets(sets, subset=subset, verbose=verbose);

  keep <- (isRawSet(sets) & isSubsetSet(sets, subset=subset));
  rawSets <- sets[keep];
  for (kk in seq(along=rawSets)) {
    verbose && enter(verbose, sprintf("Smoothing data set #%d", kk));
    set <- rawSets[[kk]];
    verbose && cat(verbose, "Name: ", set$name);

    for (h in hs) {
      hStr <- sprintf("h=%.3f", h);
      verbose && enter(verbose, hStr);
      key <- sprintf("%s,%s", set$name, hStr);
      if (key %in% names(sets)) {
        verbose && exit(verbose);
        next;
      }
      sets[[key]] <- getSmoothDataSet(set, h=h, ..., verbose=verbose);
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  sets;
} # addSmoothDataSets()


############################################################################
# HISTORY:
# 2008-07-25
# o Added file caching to getSmoothDataSet().  Great, now things are much
#   faster.
############################################################################
