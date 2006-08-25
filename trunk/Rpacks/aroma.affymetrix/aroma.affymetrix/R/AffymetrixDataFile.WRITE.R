setMethodS3("writeApd", "AffymetrixDataFile", function(this, path=NULL, ext="apd", readMap=getReadMap(this), apdMap="byChipType", ..., skip=TRUE, overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'path':
  if (is.null(path)) {
    path <- dirname(getPathname(this));
  }
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'apdMap':
  if (is.null(apdMap)) {
  } else if (identical(apdMap, "byChipType")) {
    chipType <- getChipType(this);

    # Create an optimal read map
    apdMap <- ApdMap$fromCdf(chipType=chipType);

    # Write map to current directory if non existing
    if (is.null(findApdMap(apdMap))) {
      write(apdMap);
    }
  } else if (!inherits(apdMap, "ApdMap")) {
    throw("Argument 'apdMap' is not an ApdMap object: ", class(apdMap)[1]);
  }

  verbose && enter(verbose, "Writing APD file from array '", 
                                                         getName(this), "'");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Creating APD file name
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating APD pathname");
  apdName <- paste(getName(this), ext, sep=".");
  apdName <- Arguments$getWritablePathname(apdName, path=path, 
                                         mustNotExist=(!overwrite && !skip));
  verbose && cat(verbose, "Pathname: ", apdName);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already exists?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isFile(apdName)) {
    if (skip) {
      res <- AffymetrixApdFile(apdName);
      setApdMap(res, apdMap);
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the optimal write map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(apdMap)) {
    writeMap <- getWriteMap(apdMap);
    mapType <- getMapType(apdMap);
  } else {
    writeMap <- NULL;
    mapType <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading intensities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfCells <- nbrOfProbes(this);
  chipType <- getChipType(this);
  verbose && enter(verbose, "Getting intensities from data file");
  # This may or may nottransform signals!

  y <- getProbeIntensities(this, readMap=readMap, verbose=verbose);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing APD file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating APD file");
  writeApd(apdName, data=y, chipType=chipType, mapType=mapType, writeMap=writeMap, ...);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Returning a AffymetrixApdFile object (for conveniency)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- AffymetrixApdFile(apdName);
  setApdMap(res, apdMap);

  verbose && exit(verbose);

  invisible(res);
})




############################################################################
# HISTORY:
# 2006-05-15
# o Extracted from AffymetrixDataFile.R.
# 2006-03-03
# o Added writeApd().  For now it can only write 'intensities'.
############################################################################
