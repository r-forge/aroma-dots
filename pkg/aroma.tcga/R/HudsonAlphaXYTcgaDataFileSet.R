setConstructorS3("HudsonAlphaXYTcgaDataFileSet", function(...) {
  extend(TcgaDataFileSet(...), "HudsonAlphaXYTcgaDataFileSet");
})



setMethodS3("exportTotalAndFracB", "HudsonAlphaXYTcgaDataFileSet", function(this, tags=c("*", "XY"), unf, ..., rootPath="totalAndFracBData", verbose=FALSE) {
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

    dsKK <- exportTotalAndFracB(df, dataSet=dataSet, unf=unf, 
                           rootPath=rootPath, ..., verbose=less(verbose, 5));
    verbose && cat(verbose, "Written data:");
    verbose && print(verbose, dsKK);

    verbose && exit(verbose);
  } # for (cc ...)

  # Validating everything
  res <- list();
  res$total <- AromaUnitTotalCnBinarySet$byPath(path);
  res$fracB <- AromaUnitFracBCnBinarySet$byPath(path);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2009-11-06
# o Created.
############################################################################
