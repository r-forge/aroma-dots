setConstructorS3("AromaUflFile", function(...) {
  extend(AromaAnnotationFile(...), "AromaUflFile",
    .cdf = NULL
  );
})


setMethodS3("getChipType", "AromaUflFile", function(this, ...) {
  getName(this, ...);
})

setMethodS3("getCdf", "AromaUflFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})

setMethodS3("nbrOfUnits", "AromaUflFile", function(this, ...) {
  nrow(this, ...);
})

setMethodS3("findByChipType", "AromaUflFile", function(static, chipType, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Exclude all chip type tags
  chipType <- gsub(",.*", "", chipType);

  pattern <- paste("^", chipType, "(,[.]*)*", "[.](ufl|UFL)$", sep="");
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- paste("^", chipType, "[.](ufl|UFL)[.]lnk$", sep="");
    args <- list(chipType=chipType, ...);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("fromChipType", "AromaUflFile", function(static, chipType, ...) {
  # Locate UGP file
  pathname <- findByChipType(static, chipType=chipType, ...);
  if (is.null(pathname)) {
    throw("Could not locate UFL file for this chip type: ", chipType);
  }

  # Create object
  AromaUflFile(pathname);
}, static=TRUE)



setMethodS3("readData", "AromaUflFile", function(this, ...) {
  data <- NextMethod("readData", this, ...);
  # Interpret zeros as NAs
  if (ncol(data) > 0) {
    nas <- (data[,1] == 0);
    data[nas,1] <- NA;
  }
  data;
})


setMethodS3("range", "AromaUflFile", function(this, ...) {
  lapply(this, FUN=range, ...)[[1]];
})


setMethodS3("createFromCdf", "AromaUflFile", function(static, cdf, path=getPath(cdf), tags=NULL, ...) {
  chipType <- getChipType(cdf);
  fullname <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.ufl", fullname);
  create(static, filename=filename, path=path, nbrOfRows=nbrOfUnits(cdf), types="integer", sizes=2, ...);
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
# 2007-09-11
# o Created.
############################################################################
