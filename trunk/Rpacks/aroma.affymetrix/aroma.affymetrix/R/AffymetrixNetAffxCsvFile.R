setConstructorS3("AffymetrixNetAffxCsvFile", function(..., .verify=TRUE) {
  this <- extend(AffymetrixCsvFile(..., .verify=FALSE), "AffymetrixNetAffxCsvFile");

  if (.verify)
    verify(this, ...);
  this;
})

setMethodS3("findByChipType", "AffymetrixNetAffxCsvFile", function(static, chipType, tags=".*", pattern=sprintf("^%s%s([.]|_)annot[.](csv|CSV)$", chipType, tags), ...) {
  findByChipType.AffymetrixCsvFile(static, chipType=chipType, pattern=pattern, ...);
}, static=TRUE, protected=TRUE)



setMethodS3("readDataUnitChromosomePosition", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probe[sS]etI[dD]$"="character", "^chromosome$"="character", "^(physicalPosition|chromosomeStart|probeStartPosition)$"="character"), con=NULL, ..., verbose=FALSE) {
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

  data <- readData(this, colClassPatterns=colClassPatterns, ..., verbose=less(verbose));

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



setMethodS3("readDataUnitFragmentLength", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probe[sS]etI[dD]$"="character", "^fragment.*Length.*"="character"), enzymes=1, con=NULL, ..., verbose=FALSE) {
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

  data <- readData(this, colClassPatterns=colClassPatterns, ..., verbose=less(verbose));

  # Extract fragment lengths
  verbose && enter(verbose, "Extracting fragment lengths from ([enzyme], lengths, start, stop)");

  cc <- grep("^fragment.*Length", colnames(data))[1];
  fln <- data[[cc]];
  data[[cc]] <- NULL;
  # Remove all white spaces
  fln <- gsub(" ", "", fln, fixed=TRUE);
  # Replace '///' with ';'
  fln <- gsub("///", ";", fln, fixed=TRUE);
  # Replace '//' with ','
  fln <- gsub("//", ",", fln, fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  # Are enzyme names specified?
  verbose && enter(verbose, "Inferring if enzyme names are specified");
  hasNames <- NA;
  for (kk in seq(along=fln)) {
    unit <- fln[kk]; 
    if (nchar(unit) > 0) {
      hasNames <- (regexpr("^[a-zA-Z]", unit) != -1);
      break;
    }
  }
  if (is.na(hasNames))
    throw("INTERNAL ERROR: Failed to parse CSV file for fragment lengths.");
  verbose && cat(verbose, "Has enzyme names: ", hasNames);
  verbose && exit(verbose);
  
  # Split by enzymes first
  fln <- strsplit(fln, split=";", fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  if (hasNames) {
    verbose && enter(verbose, "Identifying number of enzymes");
    nbrOfEnzymes <- max(sapply(fln, length));
    verbose && cat(verbose, "Max number of enzymes: ", nbrOfEnzymes);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Splitting into subunits and padding with NAs");
    keep <- 1:nbrOfEnzymes;
    fln <- base::lapply(fln, FUN=function(unit) {
       # Pad missing enzymes with trailing NAs
       unit <- .subset(unit, keep);
       # Split values
       strsplit(unit, split=",", fixed=TRUE);
    });
    if (isVisible(verbose, level=-50))
      verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
    verbose && exit(verbose);
    
    # Extract the name for each enzyme
    verbose && enter(verbose, "Extracting enzyme names");
    enzymeIdxs <- base::lapply(fln, FUN=function(unit) {
      base::sapply(unit, FUN=.subset, 1, USE.NAMES=FALSE);
    });
    enzymeIdxs <- unlist(enzymeIdxs, use.names=FALSE);
    # Replace '---' with NAs
    enzymeIdxs[enzymeIdxs %in% c("---")] <- NA;
    allNames <- na.omit(sort(unique(enzymeIdxs)));
    verbose && cat(verbose, "Identified enzymes: ",
                                    paste(allNames, collapse=", "));
    # Map names to indices
    enzymeIdxs <- match(enzymeIdxs,allNames);
    # Put into an ExJ matrix
    enzymeIdxs <- matrix(enzymeIdxs, nrow=nbrOfEnzymes);
    verbose && exit(verbose);
    if (isVisible(verbose, level=-50))
      verbose && str(verbose, enzymeIdxs, level=-50);
    verbose && exit(verbose);

    # Extract the fragment length for each enzyme
    verbose && enter(verbose, "Extracting fragment lengths");
    fln <- base::lapply(fln, FUN=function(unit) {
      base::sapply(unit, FUN=.subset, 2, USE.NAMES=FALSE);
    });
    fln <- unlist(fln, use.names=FALSE);
    fln <- as.integer(fln);
    # Put into an ExJ matrix
    fln <- matrix(fln, nrow=nbrOfEnzymes);
    verbose && exit(verbose);

    verbose && enter(verbose, "Sorting data by enzyme");
    # Reorganize as an JxE matrix (transposed compared with 'fln'!)
    fln2 <- matrix(NA, nrow=ncol(fln), ncol=nbrOfEnzymes);
    for (ee in 1:nbrOfEnzymes) {
      for (rr in 1:nbrOfEnzymes) {
        # Identify all indices that have enzyme 'ee' in row 'rr'
        idxs <- which(enzymeIdxs[rr,] == ee);
        if (length(idxs) > 0)
          fln2[idxs,ee] <- fln[rr,idxs];
      }
    }
    # Keep only requested enzymes
    fln <- fln2[,enzymes];
    rm(enzymeIdxs, fln2);
    verbose && exit(verbose);
  } else {
    # Extract the fragment length for each enzyme
    verbose && enter(verbose, "Extracting fragment lengths");
    fln <- base::lapply(fln, FUN=function(unit) {
      # Keep only requested enzymes
      unit <- .subset(unit, enzymes);
      parts <- strsplit(unit, split=",", fixed=TRUE);
      base::sapply(parts, FUN=.subset, 1, USE.NAMES=FALSE);
    });
    verbose && exit(verbose);
    
    # Reorganize as an JxE integer matrix
    fln <- unlist(fln, use.names=FALSE);
    fln <- as.integer(fln);
    fln <- matrix(fln, ncol=nbrOfEnzymes, byrow=TRUE);
  }
	
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
    names[-1] <- sprintf("%s.%02d", names[-1], 2:nbrOfEnzymes);

  data <- data.frame(unitName=data[[1]], fln);
  colnames(data) <- c("unitName", names);
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2007-11-28
# o Updated readDataUnitFragmentLength() of AffymetrixNetAffxCsvFile to
#   also handle "enzyme data" columns that contain named or non-named
#   multiple enzymes.
# 2007-09-16
# o Now readDataUnitChromosomePosition() of AffymetrixNetAffxCsvFile also
#   recognizes column name 'probeStartPosition'.
# 2007-09-14
# o Added support to read fragment lengths for each enzyme.
# 2007-09-11
# o Added readDataUnitFragmentLength().
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
