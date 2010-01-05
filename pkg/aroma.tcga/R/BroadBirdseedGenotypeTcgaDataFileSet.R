setConstructorS3("BroadBirdseedGenotypeTcgaDataFileSet", function(...) {
  extend(TcgaDataFileSet(...), "BroadBirdseedGenotypeTcgaDataFileSet");
})


setMethodS3("exportGenotypeCallsAndConfidenceScores", "BroadBirdseedGenotypeTcgaDataFileSet", function(this, tags=c("*", "Birdseed"), unf, ..., rootPath="callData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- unlist(strsplit(tags, split=","), use.names=FALSE);
    tags <- tags[nchar(tags) > 0];
    idxs <- which(tags == "*");
    tags[idxs] <- getTags(this, collapse=",");
    tags <- unlist(strsplit(tags, split=","), use.names=FALSE);
  }

  # Argument 'unf':
  unf <- Arguments$getInstanceOf(unf, "UnitNamesFile");

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting ", class(this)[1]);
  dataSet <- paste(c(getName(this), tags), collapse=",");
  verbose && cat(verbose, "Full data set name: ", dataSet);

  chipType <- getChipType(unf, fullname=FALSE);
  
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Exporting to path: ", path);

  for (kk in seq(this)) {
    df <- getFile(this, kk);
    
    verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                           kk, getName(df), length(this)));

    dsKK <- exportGenotypeCallsAndConfidenceScores(df, 
                dataSet=dataSet, unf=unf, rootPath=rootPath, ..., 
                verbose=less(verbose, 5));
    verbose && cat(verbose, "Written data:");
    verbose && print(verbose, dsKK);

    verbose && exit(verbose);
  } # for (cc ...)

  # Validating everything
  res <- list();
  res$acs <- AromaUnitGenotypeCallSet$byPath(path);
  pattern <- ".*,confidenceScores.asb$";
  res$ass <- AromaUnitSignalBinarySet$byPath(path, pattern=pattern);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2009-11-02
# o Added exportGenotypeCallsAndConfidenceScores() to 
#   BroadBirdseedGenotypeTcgaDataFileSet.
# 2009-10-25
# o Created.
############################################################################
