setConstructorS3("AffymetrixTsvFile", function(...) {
  this <- extend(AffymetrixFile(...), "AffymetrixTsvFile",
    "cached:.cdf" = NULL,
    "cached:.data" = NULL
  );
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})

setMethodS3("clearCache", "AffymetrixTsvFile", function(this, ...) {
  for (ff in c(".cdf", ".data")) {
    this[[ff]] <- NULL;
  }
})

setMethodS3("getChipType", "AffymetrixTsvFile", function(this, ...) {
  getName(this);
})

setMethodS3("getCdf", "AffymetrixTsvFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})

setMethodS3("findByChipType", "AffymetrixTsvFile", function(static, chipType, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- sprintf("^%s[.]tsv$", chipType);
  pathname <- findAnnotationDataByChipType(chipType, pattern);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("fromChipType", "AffymetrixTsvFile", function(static, chipType, ...) {
  # Search for the genome information file
  pathname <- static$findByChipType(chipType, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix TSV file: ", chipType);
  newInstance(static, pathname);
})

setMethodS3("verify", "AffymetrixTsvFile", function(this, ...) {
  tryCatch({
    df <- readData(this, nrows=10);
  }, error = function(ex) {
    throw("File format error of the Affymetrix TSV file: ", 
                                                  getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)

setMethodS3("getData", "AffymetrixTsvFile", function(this, force=FALSE, ...) {
  data <- this$.data;
  if (force || is.null(data)) {
    data <- readData(this, ...);
    this$.data <- data;
  }
  data;
})

setMethodS3("readData", "AffymetrixTsvFile", function(this, ..., verbose=FALSE) {
  pathname <- getPathname(this);

  colClasses <- c(
    "probeset_id"="character",
    "chr"="character",
    "snp_pos"="integer",
    "len"="integer",
    "GC"="double",
    "gc_count"="integer"
  );

  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", ...);

  names <- colnames(df);
  names <- gsub("probeset_id", "unit", names);
  names <- gsub("chr", "chromosome", names);
  names <- gsub("GC", "gc", names);
  names <- gsub("snp_pos", "physical position", names);
  names <- gsub("_", " ", names);
  names <- toCamelCase(names);
  colnames(df) <- names;

  # Rescale GC contents to [0,1]
  df[["gc"]] <- df[["gc"]]/100;

  # Remap chromsome X->23, Y->24
  chr <- df[["chromosome"]];
  chr[chr == "X"] <- 23;
  chr[chr == "Y"] <- 24;
  suppressWarnings({
    chr <- as.integer(chr);
  })
  df[["chromosome"]] <- chr;
  rm(chr);

  # Remove duplicated rows (type by Affymetrix?!? /HB 2007-04-02)
  df <- unique(df);

  gc();

  # Convert unit names to unit indices
  cdf <- getCdf(this);
  units <- match(df[["unit"]], getUnitNames(cdf));
  if (any(is.na(units))) {
    throw("File format error: Identified units that do not exist in the CDF: ", getChipType(cdf));
  }
  df[["unit"]] <- units;

#  rownames(df) <-  units;
  
  o <- order(units);
  df <- df[o,];

  df;
})

setMethodS3("getField", "AffymetrixTsvFile", function(this, units=NULL, field, ...) {
  if (is.null(units)) {
    cdf <- getCdf(this);
    units <- 1:nbrOfUnits(cdf);
  }

  data <- getData(this, ...);
  if (!field %in% colnames(data))
    throw("No such field: ", field);

  idxs <- match(units, data$unit);
  data[[field]][idxs];
})

setMethodS3("getPosition", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="physicalPosition", ...);
})

setMethodS3("getFragmentLengths", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="len", ...);
})

setMethodS3("getGc", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="gc", ...);
})



############################################################################
# HISTORY:
# 2007-03-02
# o Created.
############################################################################
