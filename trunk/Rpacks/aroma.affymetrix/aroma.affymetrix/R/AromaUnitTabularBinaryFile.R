setConstructorS3("AromaUnitTabularBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), "AromaUnitTabularBinaryFile",
    .cdf = NULL
  );
})


setMethodS3("getFilenameExtension", "AromaUnitTabularBinaryFile", static=TRUE, abstract=TRUE);

setMethodS3("nbrOfUnits", "AromaUnitTabularBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("getChipType", "AromaUnitTabularBinaryFile", function(this, ...) {
  getName(this, ...);
})


setMethodS3("getCdf", "AromaUnitTabularBinaryFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("fromChipType", "AromaUnitTabularBinaryFile", function(static, chipType, ...) {
  pathname <- findByChipType(static, chipType=chipType, ...);
  if (is.null(pathname)) {
    throw("Could not locate file for this chip type: ", chipType);
  }

  # Create object
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("findByChipType", "AromaUnitTabularBinaryFile", function(static, chipType, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Exclude all chip type tags
  chipType <- gsub(",.*", "", chipType);

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s(,[.]*)*[.]%s$", chipType, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s(,[.]*)*[.]%s[.]lnk$", chipType, ext);
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


setMethodS3("indexOfUnits", "AromaUnitTabularBinaryFile", function(this, names, ...) {
  # Look up unit names from CDF
  cdf <- getCdf(this);
  idxs <- match(names, getUnitNames(cdf));
  idxs;
}, protected=TRUE)



setMethodS3("allocateFromCdf", "AromaUnitTabularBinaryFile", function(static, cdf, path=getPath(cdf), tags=NULL, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(cdf);
  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create tabular binary file
  allocate(static, filename=filename, path=path, nbrOfRows=nbrOfUnits(cdf), ...);
}, static=TRUE)



setMethodS3("importFrom", "AromaUnitTabularBinaryFile", function(this, object, ...) {
  if (inherits(object, "AffymetrixNetAffxCsvFile")) {
    importFromAffymetrixNetAffxCsvFile(this, object, ...);
  } else if (inherits(object, "DChipGenomeInformation")) {
    importFromDChipGenomeInformation(this, object, ...);
  } else if (inherits(object, "GenomeInformation")) {
    importFromGenomeInformation(this, object, ...);
  } else if (inherits(object, "AffymetrixTabularFile")) {
    importFromAffymetrixTabularFile(this, object, ...);
  } else {
    throw("Do not know how to import from an object of class ", 
                                                          class(object)[1]);
  }
})

setMethodS3("importFromAffymetrixTabularFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromDChipGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
