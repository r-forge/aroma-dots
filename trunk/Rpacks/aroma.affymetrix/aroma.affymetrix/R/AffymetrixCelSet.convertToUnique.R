#*/###########################################################################
setMethodS3("convertToUnique", "AffymetrixCelSet", function(this, ..., tags="UNQ", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Converting to unique CDF");
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already unique?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this)
  if (isUniqueCdf(cdf)) {
    verbose && cat(verbose, "Already based on a unique CDF");
    verbose && exit(verbose);
    return(invisible(this));
  } else {
    verbose && enter(verbose, "Getting unique CDF");
    cdfUnique <- getUniqueCdf(cdf)
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # getting output directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tagsInput <- paste(getTags(this), collapse=",")
  tags <- paste(tagsInput, tags, sep=",")
  verbose && cat(verbose, "Tags:", tags);
  outputPath <- gsub( tagsInput, tags, getPath(this) )

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # check if already done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Test whether dataset exists");
  # HB: I consider try() obsolete.
  # HB: Don't think argument 'chipType' makes a difference if 'cdf' is given.
  tryCatch({
    outputDataSet <- AffymetrixCelSet$byName(getName(this), tags=tags, 
                                            verbose=verbose, cdf=cdfUnique, 
                          chipType=getChipType(this), checkChipType=FALSE);
  }, error = function(ex) {});
  
  if (inherits(outputDataSet, "AffymetrixCelSet")) {
    verbose && cat(verbose, "Dataset already created.");
    verbose && exit(verbose);
    return(invisible(outputDataSet));
  }
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read indices for old and new
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices from standard CDF");
  cdfStandard <- readCdf(getPathname(cdf), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Reading cell indices list from unique CDF");
  cdfUniqueIndices <- readCdf(getPathname(cdfUnique), units=NULL, readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE,readUnitType=FALSE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=TRUE, readIsPm=FALSE);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(this);
  df <- getFile(this, 1);
  
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Do the conversion from standard CDF to unique CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {

      verbose && enter(verbose, "Converting CEL data from standard to unique CDF");
  
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Read data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading values according to standard CDF");
      df <- getFile(this, kk);
      data <- readCelUnits(getPathname(df), cdf=cdfStandard);
      hdr <- readCelHeader(getPathname(df));
      verbose && exit(verbose);

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Write data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Create CEL file to store results, if missing
      nr <- nbrOfRows(cdfUnique);
      nc <- nbrOfColumns(cdfUnique);
      header <- list(filename=pathname, version=4, rows=nr, cols=nc, total=nr*nc, noutliers=0, chiptype=getChipType(cdfUnique), header=hdr$header, algorithm=hdr$algorithm, parameters=hdr$parameters, cellmargin=hdr$cellmargin, nmasked=hdr$nmasked);

      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createCel(pathname, header=header);
      verbose && cat(verbose, "Writing values according to unique CDF");
      updateCelUnits(pathname, cdf=cdfUniqueIndices, data=data, verbose=FALSE);
      verbose && exit(verbose);

      rm(data,hdr);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
  } # for (kk ...)

  outputDataSet <- AffymetrixCelSet$fromName(getName(this), tags=tags, cdf=cdfUnique, verbose=verbose, checkChipType=FALSE);
  verbose && exit(verbose);
  
  invisible(outputDataSet);
})
