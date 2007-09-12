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


setMethodS3("createFromCdf", "AromaUgpFile", function(static, cdf, path=getPath(cdf), ...) {
  chipType <- getChipType(cdf);
  create(static, chipType=chipType, nbrOfElements=nbrOfUnits(cdf), path=path, ...);
}, static=TRUE)


setMethodS3("createFromGenomeInformation", "AromaUgpFile", function(static, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  chipType <- getChipType(gi);
  cdf <- AffymetrixCdfFile$fromChipType(chipType);
  ugp <- createFromCdf(cdf, ...);
  importFromGenomeInformation(ugp, gi);
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
    chr[chr == "Y"] <- 24;
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



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUgpFile", function(this, csv, shift=0, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  if (!inherits(csv, "AffymetrixNetAffxCsvFile")) {
    throw("Argument 'csv' is not an AffymetrixNetAffxCsvFile: ", class(csv)[1]);
  }

  # Argument 'shift':
  shift <- Arguments$getInteger(shift);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing (unit name, chromosome, position) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);
  cdfUnitNames <- getUnitNames(cdf);

  # Read data
  data <- readDataUnitChromosomePosition(csv, ..., verbose=less(verbose));

  # Map to CDF unit names
  cdfUnits <- match(data[[1]], cdfUnitNames);

  # Exclude units that are not in the CDF
  keep <- which(!is.na(cdfUnits));
  cdfUnits <- cdfUnits[keep];
  if (length(cdfUnits) == 0) {
    warning("None of the imported unit names match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }
  data <- data[keep,2:3,drop=FALSE];
  if (shift != 0) {
    data[[2]] <- data[[2]] + shift;
  }
 
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  res <- updateData(this, idxs=cdfUnits, chromosome=data[[1]], position=data[[2]], verbose=less(verbose));

  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(cdfUnits);
})


############################################################################
# HISTORY:
# 2007-09-10
# o Added importFromAffymetrixNetAffxCsvFile() to AromaUgpFile.
# 2007-03-04
# o Added findByChipType() and fromChipType().
# o Now the default path for createFromCdf() is the same as for the CDF.
# 2007-03-03
# o Now inherits from generic AromaGenomePositionFile.
# 2007-03-02
# o Created. Can import genome information data.
############################################################################
