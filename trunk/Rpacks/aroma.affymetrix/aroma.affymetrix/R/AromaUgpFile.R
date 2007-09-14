setConstructorS3("AromaUgpFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUgpFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("getFilenameExtension", "AromaUgpFile", function(static, ...) {
  "ugp";
}, static=TRUE)


setMethodS3("readData", "AromaUgpFile", function(this, ...) {
  data <- NextMethod("readData", this, ...);

  # Interpret zeros as NAs
  if (ncol(data) > 0) {
    for (cc in seq(length=ncol(data))) {
      nas <- (data[,cc] == 0);
      data[nas,cc] <- NA;
    }
  }

  data;
})

setMethodS3("getGenomeVersion", "AromaUgpFile", function(this, ...) {
  tags <- getTags(this, ...);
  tags <- grep("^hg", tags, value=TRUE);
  tags;
}, protected=TRUE)


setMethodS3("allocateFromCdf", "AromaUgpFile", function(static, ...) {
  # NextMethod() not supported here.
  allocateFromCdf.AromaUnitTabularBinaryFile(static, ..., types=rep("integer",2), sizes=c(1,4));
}, static=TRUE)



setMethodS3("getUnitsAt", "AromaUgpFile", function(this, chromosome, range=NULL, ..., verbose=FALSE) {
  # Stratify by chromosome
  data <- this[,1,drop=TRUE];
  keep <- !is.na(data) & (data %in% chromosome);
  idxs <- which(keep);

  if (!is.null(range)) {
    data <- this[idxs,2,drop=TRUE];
    keep <- !is.na(data);
    keep <- keep & (range[1] <= data & data <= range[2]);
    idxs <- idxs[keep];
  }
  
  idxs;
}, protected=TRUE)





setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUgpFile", function(this, csv, shift="auto", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  if (!inherits(csv, "AffymetrixNetAffxCsvFile")) {
    throw("Argument 'csv' is not an AffymetrixNetAffxCsvFile: ", class(csv)[1]);
  }

  # Argument 'shift':
  if (!identical(shift, "auto"))
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
  importNames <- attr(data, "importNames");

  # Map to CDF unit names
  cdfUnits <- match(data[,1], cdfUnitNames);

  # Exclude units that are not in the CDF
  keep <- which(!is.na(cdfUnits));
  cdfUnits <- cdfUnits[keep];
  if (length(cdfUnits) == 0) {
    warning("None of the imported unit names match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }
  data <- data[keep,2:3,drop=FALSE];
  importNames <- importNames[2:3];

  # Shift positions?
  if (identical(shift, "auto")) {
    shift <- 0;
    if ("chromosomeStart" %in% importNames)
      shift <- 13;
  }
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.", shift);
    data[,2] <- data[,2] + as.integer(shift);
  }
 
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  this[cdfUnits,1] <- data[,1];
  this[cdfUnits,2] <- data[,2];

  rm(data);
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(cdfUnits);
})



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
    chr[chr == "Z"] <- 25;
    suppressWarnings({
      chr <- as.integer(chr);
    })
  }
  
  pos <- data[,"physicalPosition"];
  suppressWarnings({
    pos <- as.integer(pos);
  })

  this[,1] <- chr;
  this[,2] <- pos;
})



############################################################################
# HISTORY:
# 2007-09-13
# o Removed createFromGenomeInformation().
# o Updated AromaUgpFile according to changes in super class.
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
