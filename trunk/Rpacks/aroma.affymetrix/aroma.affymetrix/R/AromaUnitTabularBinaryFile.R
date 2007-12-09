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


setMethodS3("fromChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  pathname <- findByChipType(static, chipType=chipType, tags=tags, ...);
  if (is.null(pathname)) {
    throw("Could not locate file for this chip type: ", 
                                   paste(c(chipType, tags), collapse=","));
  }

  # Create object
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("findByChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  chipType <- parts[1];
  tags <- parts[-1];

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s.*[.]%s[.]lnk$", chipType, ext);
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



setMethodS3("importFrom", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  if (inherits(src, "AffymetrixNetAffxCsvFile")) {
    importFromAffymetrixNetAffxCsvFile(this, src, ...);
  } else if (inherits(src, "DChipGenomeInformation")) {
    importFromDChipGenomeInformation(this, src, ...);
  } else if (inherits(src, "GenomeInformation")) {
    importFromGenomeInformation(this, src, ...);
  } else if (inherits(src, "AffymetrixTabularFile")) {
    importFromAffymetrixTabularFile(this, src, ...);
  } else if (inherits(src, "GenericTabularFile")) {
    importFromGenericTabularFile(this, src, ...);
  } else {
    throw("Do not know how to import from an src of class ", class(src)[1]);
  }
})


setMethodS3("importFromGenericTabularFile", "AromaUnitTabularBinaryFile", abstract=TRUE);

setMethodS3("importFromAffymetrixTabularFile", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  if (!inherits(src, "AffymetrixTabularFile")) {
    throw("Argument 'src' is not a AffymetrixTabularFile file: ", class(src)[1]);
  }

  importFromGenomeInformation(this, src, ...);
});

setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromDChipGenomeInformation", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  if (!inherits(src, "DChipGenomeInformation")) {
    throw("Argument 'src' is not a DChipGenomeInformation file: ", class(src)[1]);
  }

  importFromGenomeInformation(this, src, ...);
})


setMethodS3("importFromGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
