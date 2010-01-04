setConstructorS3("HudsonAlphaGenotypeTcgaDataFileSet", function(...) {
  extend(TcgaDataFileSet(...), "HudsonAlphaGenotypeTcgaDataFileSet");
})


setMethodS3("exportGenotypeCalls", "HudsonAlphaGenotypeTcgaDataFileSet", function(this, tags=c("*", "BeadStudio"), unf, ..., rootPath="callData", verbose=FALSE) {
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

  for (kk in seq(this)) {
    df <- getFile(this, kk);
    
    verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                           kk, getName(df), length(this)));

    dsKK <- exportGenotypeCalls(df, 
                dataSet=dataSet, unf=unf, rootPath=rootPath, ..., 
                verbose=less(verbose, 5));
    verbose && cat(verbose, "Written data:");
    verbose && print(verbose, dsKK);

    verbose && exit(verbose);
  } # for (cc ...)

  # Validating everything
  res <- list();
  res$acs <- AromaUnitGenotypeCallSet$byPath(path);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2009-12-05
# o Created.
############################################################################
