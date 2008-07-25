# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Public functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
addRocDataSet <- function(sets, key, name=key, ...) {
  # Already added?
  if (name %in% names(sets))
    return(sets);

  cat("Adding data set: ", name, "\n", sep="");

  # Assert that it exists
  pathname <- getRocPathname(key);
  set <- list(name=name, pathname=pathname);
  sets[[name]] <- set;

  sets;
}



addRocData <- function(set, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  res <- set$rocData;
  if (!force && !is.null(res))
    return(set);

  verbose && enter(verbose, "Building RocData");
  verbose && cat(verbose, "Data set: ", set$name);

  # Get log-ratios
  M <- loadLogRatios(set);

  # Get true CNs
  C <- set$truth;
  if (is.null(C)) {
    verbose && enter(verbose, "Creating default truth (by columns)");
    C <- rep(1, times=ncol(M));
    C[(n23[colnames(M)] == 2)] <- 0;
    C <- as.integer(C);
    verbose && str(verbose, C);
    verbose && exit(verbose);
  }

  res <- RocData(truth=C, data=M);
  verbose && exit(verbose);

  set$rocData <- res;

  set;
} # addRocData()



removeRawDataSets <- function(sets, ...) {
  sets[!isRawDataSet(sets)];
}
