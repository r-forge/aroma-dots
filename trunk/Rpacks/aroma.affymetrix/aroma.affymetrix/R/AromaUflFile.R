setConstructorS3("AromaUflFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUflFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("getFilenameExtension", "AromaUflFile", function(static, ...) {
  "ufl";
}, static=TRUE, protected=TRUE);


setMethodS3("nbrOfEnzymes", "AromaUflFile", function(this, ...) {
  nbrOfColumns(this, ...);
})


setMethodS3("readData", "AromaUflFile", function(this, ...) {
  data <- NextMethod("readData", this, ...);

  # Interpret zeros as NAs
  if (ncol(data) > 0) {
    nas <- (data[,1] == 0);
    data[nas,1] <- NA;
  }

  data;
})

setMethodS3("allocateFromCdf", "AromaUflFile", function(static, cdf, nbrOfEnzymes=1, ...) {
  # Argument 'nbrOfEnzymes':
  nbrOfEnzymes <- Arguments$getInteger(nbrOfEnzymes, range=c(1,10));

  types <- rep("integer", nbrOfEnzymes);
  sizes <- rep(2, nbrOfEnzymes);

  # NextMethod() not supported here.
  allocateFromCdf.AromaUnitTabularBinaryFile(static, cdf=cdf, types=types, sizes=sizes, ...);
}, static=TRUE)



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUflFile", function(this, csv, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  if (!inherits(csv, "AffymetrixNetAffxCsvFile")) {
    throw("Argument 'csv' is not an AffymetrixNetAffxCsvFile: ", class(csv)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing (unit name, fragment length) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);
  cdfUnitNames <- getUnitNames(cdf);

  # Read data
  data <- readDataUnitFragmentLength(csv, ..., verbose=less(verbose));

  # Map to CDF unit names
  cdfUnits <- match(data[[1]], cdfUnitNames);

  # Exclude units that are not in the CDF
  keep <- which(!is.na(cdfUnits));
  cdfUnits <- cdfUnits[keep];
  if (length(cdfUnits) == 0) {
    warning("None of the imported unit names match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }
  data <- data[keep,2,drop=TRUE];
 
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  # Update
  this[cdfUnits,1] <- data;
  rm(data);

  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(cdfUnits);
})


############################################################################
# HISTORY:
# 2007-09-14
# o Added support for multiple fragment lengths, in case multiple enzymes
#   were used for the same assay.
# 2007-09-11
# o Created.
############################################################################
