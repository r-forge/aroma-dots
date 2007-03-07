setConstructorS3("AromaUgpFile", function(...) {
  this <- extend(AromaGenomePositionFile(...), "AromaUgpFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("getFilenameExtension", "AromaUgpFile", function(static, ...) {
  "ugp";
}, static=TRUE)

setMethodS3("nbrOfUnits", "AromaUgpFile", function(this, ...) {
  nbrOfElements(this, ...);
})


setMethodS3("findByChipType", "AromaUgpFile", function(static, chipType, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Exclude all chip type tags
  chipType <- gsub(",.*", "", chipType);

  pattern <- paste("^", chipType, "(,[.]*)*", "[.](ugp|UGP)$", sep="");
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- paste("^", chipType, "[.](c|C)(d|D)(f|F)[.]lnk$", sep="");
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


setMethodS3("fromChipType", "AromaUgpFile", function(static, chipType, ...) {
  # Locate UGP file
  pathname <- findByChipType(static, chipType=chipType, ...);
  if (is.null(pathname)) {
    throw("Could not locate UGP file for this chip type: ", chipType);
  }

  # Create object
  AromaUgpFile(pathname);
}, static=TRUE)

setMethodS3("fromCdf", "AromaUgpFile", function(static, cdf, ...) {
  fromChipType(static, getChipType(cdf));
})

setMethodS3("getCdf", "AromaUgpFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  
  cdf;
})

setMethodS3("createFromCdf", "AromaUgpFile", function(static, cdf, path=getPath(cdf), ...) {
  chipType <- getChipType(cdf);
  create(static, chipType=chipType, nbrOfElements=nbrOfUnits(cdf), path=path, ...);
}, static=TRUE)


setMethodS3("indexOfElements", "AromaUgpFile", function(this, names, ...) {
  # Look up unit names from CDF
  cdf <- getCdf(this);
  idxs <- match(names, getUnitNames(cdf));
  idxs;
}, protected=TRUE)


setMethodS3("getUnitsAt", "AromaUgpFile", function(this, ...) {
  getElementsAt(this, ...);
}, protected=TRUE)


setMethodS3("importFromGenomeInformation", "AromaUgpFile", function(this, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  # AD HOC patch, since units==NULL does not work./HB 2007-03-03
  units <- seq_len(nbrOfUnits(gi));
  data <- getData(gi, units=units, fields=c("chromosome", "physicalPosition"));

  chr <- data[,"chromosome"];
  if (is.character(chr)) {
    chr[chr == "X"] <- 23;
    chr[chr == "Y"] <- 23;
    suppressWarnings({
      chr <- as.integer(chr);
    })
  }
  
  pos <- data[,"physicalPosition"];
  suppressWarnings({
    pos <- as.integer(pos);
  })

  updateData(this, chromosome=chr, position=pos);
})


setMethodS3("createFromGenomeInformation", "AromaUgpFile", function(static, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  chipType <- getChipType(gi);
  cdf <- AffymetrixCdfFile$fromChipType(chipType);
  ugp <- createFromCdf(cdf, ...);
  importFromGenomeInformation(ugp, gi);
}, static=TRUE)



############################################################################
# HISTORY:
# 2007-03-04
# o Added findByChipType() and fromChipType().
# o Now the default path for createFromCdf() is the same as for the CDF.
# 2007-03-03
# o Now inherits from generic AromaGenomePositionFile.
# 2007-03-02
# o Created. Can import genome information data.
############################################################################
