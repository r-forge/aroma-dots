setConstructorS3("AffymetrixNetAffxCsvFile", function(..., .verify=TRUE) {
  this <- extend(AffymetrixCsvFile(..., .verify=FALSE), "AffymetrixNetAffxCsvFile");

  if (.verify)
    verify(this, ...);
  this;
})

setMethodS3("findByChipType", "AffymetrixNetAffxCsvFile", function(static, chipType, tags=".*", pattern=sprintf("^%s%s([.]|_)annot[.](csv|CSV)$", chipType, tags), ...) {
  findByChipType.AffymetrixCsvFile(static, chipType=chipType, pattern=pattern, ...);
}, static=TRUE, protected=TRUE)



setMethodS3("readDataUnitChromosomePosition", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^(physicalPosition|chromosomeStart)$"="character"), con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading (unitName, fragmentLength) from file");

  data <- readData(this, colClassPatterns=colClassPatterns, camelCaseNames=TRUE, ..., verbose=less(verbose));

  # Convert chromosome strings to integers
  cc <- grep("^chr", colnames(data))[1];
  map <- c(X=23, Y=24, Z=25);
  for (kk in seq(along=map)) {
    data[[cc]] <- gsub(names(map)[kk], map[kk], data[[cc]]);
  }
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  # Convert positions to integers
  cc <- grep("(p|P)os", colnames(data));
  if (length(cc) == 0)
    cc <- grep("(s|S)tart", colnames(data));
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  colnames(data) <- c("unitName", "chromosome", "position");
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



setMethodS3("readDataUnitFragmentLength", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^fragmentLength.*"="character"), con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading (unitName, fragmentLength) from file");

  data <- readData(this, colClassPatterns=colClassPatterns, camelCaseNames=TRUE, ..., verbose=less(verbose));

  # Extract fragment lengths
  verbose && enter(verbose, "Extracting fragment lengths from (lengths, start, stop)");

  cc <- grep("^fragmentLength", colnames(data))[1];
  fln <- data[[cc]];
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln, level=-50);
  fln <- strsplit(fln, split="//", fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln, level=-50);
  fln <- sapply(fln, FUN=.subset, 1);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln, level=-50);
  fln <- trim(fln);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln, level=-50);
  fln <- as.integer(fln);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln, level=-50);
  data[[cc]] <- fln;
  verbose && exit(verbose);

  gc <- gc();

  colnames(data) <- c("unitName", "fragmentLength");
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2007-09-11
# o Added readDataUnitFragmentLength().
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
