setConstructorS3("AffymetrixCsvGenomeInformation", function(...) {
  this <- extend(GenomeInformation(...), "AffymetrixCsvGenomeInformation");
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})

setMethodS3("findByChipType", "AffymetrixCsvGenomeInformation", function(static, chipType, version=NULL, ...) {
  # Argument 'version':
  if (is.null(version))
    version <- ".*";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- sprintf("^%s[.].*[.]annot[.]csv$", chipType);
  pathname <- findAnnotationDataByChipType(chipType, pattern);

  pathname;
}, static=TRUE, protected=TRUE)

setMethodS3("fromChipType", "AffymetrixCsvGenomeInformation", function(static, chipType, version=NULL, ...) {
  # Search for the genome information file
  pathname <- static$findByChipType(chipType, version=version, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix CSV annotation file: ", chipType);
  newInstance(static, pathname);
})

setMethodS3("verify", "AffymetrixCsvGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readData(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the Affymetrix CSV annotation file: ", 
                                                  getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)

setMethodS3("readData", "AffymetrixCsvGenomeInformation", function(this, ..., verbose=FALSE) {
  pathname <- getPathname(this);
  hdr <- scan(pathname, what=character(0), nlines=1, sep=",", quote="\"", quiet=TRUE);
  nbrOfColumns <- length(hdr);
  colClasses <- rep("NULL", nbrOfColumns);
  names(colClasses) <- hdr;
  colClasses["Probe Set ID"] <- "character";
  colClasses["Chromosome"] <- "character";
  colClasses["Physical Position"] <- "character";
  df <- read.table(pathname, colClasses=colClasses, header=TRUE, sep=",", quote="\"", fill=TRUE, check.names=FALSE, na.strings=c("---"), ...);
  df[["Physical Position"]] <- as.integer(df[["Physical Position"]]);
  colnames(df) <- toCamelCase(colnames(df));
  
  df;
})


############################################################################
# HISTORY:
# 2007-03-02
# o Created.
############################################################################
