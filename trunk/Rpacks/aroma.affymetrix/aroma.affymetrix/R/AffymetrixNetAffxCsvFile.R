setConstructorS3("AffymetrixNetAffxCsvFile", function(..., .verify=TRUE) {
  this <- extend(AffymetrixCsvFile(..., .verify=FALSE), "AffymetrixNetAffxCsvFile");

  if (.verify)
    verify(this, ...);
  this;
})

setMethodS3("findByChipType", "AffymetrixNetAffxCsvFile", function(static, chipType, tags=".*", pattern=sprintf("^%s%s([.]|_)annot[.](csv|CSV)$", chipType, tags), ...) {
  findByChipType.AffymetrixCsvFile(static, chipType=chipType, pattern=pattern, ...);
}, static=TRUE, protected=TRUE)



setMethodS3("readDataUnitChromosomePosition", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^chromosome$"="character", "^(physicalPosition|chromosomeStart)$"="character"), con=NULL, ..., verbose=FALSE) {
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
  if (length(cc) == 0 || is.na(cc)) {
    throw("Failed to locate chromosome column.");
  }
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
  if (length(cc) == 0 || is.na(cc))
    throw("Failed to locate position column.");
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  attr(data, "importNames") <- colnames(data);
  colnames(data) <- c("unitName", "chromosome", "position");
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



setMethodS3("readDataUnitFragmentLength", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^fragmentLength.*"="character"), enzymes=1, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfEnzymes':
  enzymes <- Arguments$getIndices(enzymes, range=c(1,10));
  if (any(duplicated(enzymes))) {
    throw("Argument 'enzymes' contains duplicated values: ", 
                                          paste(enzymes, collapse=", "));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfEnzymes <- length(enzymes);

  verbose && enter(verbose, "Reading (unitName, fragmentLength+) from file");

  data <- readData(this, colClassPatterns=colClassPatterns, camelCaseNames=TRUE, ..., verbose=less(verbose));

  # Extract fragment lengths
  verbose && enter(verbose, "Extracting fragment lengths from (lengths, start, stop)");

  cc <- grep("^fragmentLength", colnames(data))[1];
  fln <- data[[cc]];
  data[[cc]] <- NULL;
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  # Split by enzymes first
  fln <- strsplit(fln, split="///", fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  # Keep only requested enzymes
  fln <- base::lapply(fln, FUN=.subset, enzymes);
  gc <- gc();

  # Extract the fragment length for each enzyme
  fln <- base::lapply(fln, FUN=function(unit) {
    parts <- strsplit(unit, split="//", fixed=TRUE);
    value <- base::sapply(parts, FUN=.subset, 1, USE.NAMES=FALSE);
    value <- trim(value);
    value;
  })
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  # Turn into an integer matrix
  fln <- unlist(fln, use.names=FALSE);
  fln <- as.integer(fln);
  gc <- gc();
  fln <- matrix(fln, ncol=nbrOfEnzymes, byrow=TRUE);
  gc <- gc();
	
  if (isVisible(verbose, level=-10))
    verbose && str(verbose, fln, level=-10);

  # Sanity check
  if (nrow(fln) != nrow(data)) {
    throw("Internal error. nrow(fln) != nrow(data): ", 
                                       nrow(fln), " != ", nrow(data));
  }

  verbose && exit(verbose);

  names <- c("fragmentLength", rep(c("fragmentLength"), nbrOfEnzymes-1));
  if (nbrOfEnzymes > 1)
    names[-1] <- sprintf("%s,%02d", names[-1], 2:nbrOfEnzymes);

  data <- data.frame(unitName=data[[1]], fln);
  colnames(data) <- c("unitName", names);
  attr(data, "header") <- NULL;

  gc <- gc();

  verbose && exit(verbose);

  data;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2007-09-14
# o Added support to read fragment lengths for each enzyme.
# 2007-09-11
# o Added readDataUnitFragmentLength().
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
