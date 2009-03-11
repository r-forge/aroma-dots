setMethodS3("findBatchesByName", "AffymetrixCelSet", function(static, name, tags=NULL, chipType, pattern=sprintf("^%s$", paste(c(name, "[0-9]+", tags), collapse=",")), ..., paths=c("rawData", "probeData"), verbose=FALSE) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);
  
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);
  
  # Argument ' pattern':
  pattern <- Arguments$getRegularExpression(pattern);
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Searching for a batch collect of data sets");
  verbose && cat(verbose, "Data set: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=","));
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Pattern: ", pattern);

  # List all file sets
  path <- paths[1];
  paths <- list.files(pattern=pattern, path=path, full.names=TRUE);
  # Nothing to do?
  if (length(paths) == 0) {
    return(NULL);
  }

  # Keep only directories
  paths <- paths[sapply(paths, FUN=isDirectory)];
  # Nothing to do?
  if (length(paths) == 0) {
    return(NULL);
  }
  verbose && cat(verbose, "Found possible data sets:");
  verbose && str(verbose, paths);

  # Scan for the correct chip type
  paths <- file.path(paths, chipType);
  verbose && str(verbose, paths);
  paths <- sapply(paths, FUN=function(path) {
    Arguments$getReadablePath(path);
  });
  verbose && str(verbose, paths);

  # Keep only directories
  paths <- paths[sapply(paths, FUN=isDirectory)];
  # Nothing to do?
  if (length(paths) == 0) {
    return(NULL);
  }
  verbose && str(verbose, paths);

  # Order in lexicographic order
  paths <- sort(paths);
  names(paths) <- NULL;

  verbose && cat(verbose, "Located ", length(paths), " data sets:");
  dataSets <- basename(dirname(paths));
  verbose && print(verbose, dataSets);

  verbose && exit(verbose);
  
  paths;
}, static=TRUE)





setMethodS3("batchesByName", "AffymetrixCelSet", function(static, name, tags=NULL, chipType=NULL, ..., cdf=NULL, merge=TRUE, verbose=FALSE) {
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Searching for a batches by name");
  if (!is.null(cdf)) {
    chipType <- getChipType(cdf, fullname=FALSE);
  }
  
  verbose && enter(verbose, "Searching for a sequence of batch data sets");
  pathnames <- findBatchesByName(static, name=dataSet, tags=tags, chipType=chipType, verbose=less(verbose, 10));
  verbose && str(verbose, pathnames);
  verbose && exit(verbose);

  dataSets <- dirname(pathnames);
  dataSets <- basename(dataSets);
  print(dataSets);


  if (merge) {
    res <- NULL;
  } else {
    res <- list();
  }

  for (kk in seq(along=dataSets)) {
    dataSet <- dataSets[kk];
    verbose && enter(verbose, sprintf("Batch #%d of %d", kk, length(dataSets)));
    verbose && cat(verbose, "Data set: ", dataSet);
    verbose && cat(verbose, "Tags: NULL");
    verbose && cat(verbose, "Chip type: ", chipType);

#    ds <- static$fromFiles(pathname, ..., verbose=verbose);
    ds <- byName(static, name=dataSet, tags=NULL, chipType=chipType, ..., verbose=verbose);
    print(ds);

    if (merge) {
      if (is.null(res)) {
        res <- ds;
      } else {
        res <- append(res, ds);
        setTags(res, intersect(getTags(res), getTags(ds)));
      }
    } else {
      res[[kk]] <- ds;
    }
    rm(ds);
    verbose && exit(verbose);
  } # for (kk ...)

  print(res);

  verbose && exit(verbose);

  res;
}, static=TRUE) # batchesByName()


############################################################################
# HISTORY:
# 2009-02-25 [HB]
# o Created for simple setup and merging of batch sets (as PN have).
############################################################################
