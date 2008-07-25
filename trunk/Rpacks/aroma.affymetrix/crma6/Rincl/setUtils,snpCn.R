# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Public functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getSubsetDataSet <- function(set, subset=c("snp", "cn"), ..., verbose=FALSE) {
  # Argument 'subset':
  subset <- match.arg(subset);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting subset");
  verbose && cat(verbose, "Data set: ", set$name);
  verbose && cat(verbose, "Subset: ", subset);
  name <- sprintf("%s,subset=%s", set$name, subset);
  verbose && cat(verbose, "New data set: ", name);

  # Subset by unit type?
  types <- attr(unitTypes, "map");
  keep <- whichVector(unitTypes == types[subset]);

  # Get log-ratios
  set <- addRocData(set, ...);
  rocData <- extractSubset(set$rocData, rows=keep);
  set$rocData <- rocData;
  rm(rocData);

  set$name <- name;

  # Get positions
  if (!is.null(set$positions))
    set$positions <- set$positions[keep];
  rm(keep);

  verbose && exit(verbose);

  set;
} # getSubsetDataSet()



addSubsetDataSets <- function(sets, subset, ..., verbose=FALSE) {
  # Nothing to do?
  if (subset == "all")
    return(sets);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  rawSets <- sets[isRawSet(sets) & !isSubsetSet(sets)];
  for (kk in seq(along=rawSets)) {
    set <- rawSets[[kk]];
    key <- sprintf("%s,subset=%s", set$name, subset);
    if (key %in% names(sets))
      next;

    verbose && enter(verbose, sprintf("Subsetting data set %d", kk));
    set <- getSubsetDataSet(set, subset=subset, ..., verbose=verbose);

    sets[[key]] <- set;
    rm(set);

    verbose && exit(verbose);
  }

  sets;
} # addSubsetDataSets()
