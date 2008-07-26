# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Public functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getSmoothDataSet <- function(set, h, ..., verbose=FALSE) {
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

  # Sanity check
  if (is.null(set$rocData)) {
    throw("Cannot get smoothed data set, because the original data set has yet to be loaded.");
  }

  # Clone data set object  
  res <- set;

  # Update the name
  res$name <- sprintf("%s,h=%.3f", set$name, h);

  # Smooth the ROC data
  res$rocData <- extractSmoothRocData(res$rocData, h=h, verbose=verbose);

  verbose && exit(verbose);

  res;
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
