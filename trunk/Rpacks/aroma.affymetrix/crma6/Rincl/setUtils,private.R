# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Private functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getRocPathname <- function(tags, ...) {
  # Argument 'tags':
  tags <- Arguments$getCharacter(tags);

  name <- dataSetName;
  tags <- paste(c(name, tags), collapse=",");
  filenames <- c(tags, sprintf("%s,M", tags));
  filenames <- paste(filenames, ",units23.RData", sep="");
  pathnames <- file.path(dataPath, filenames);
  pathnames <- sapply(pathnames, FUN=filePath, expandLinks="any");
  keep <- sapply(pathnames, FUN=isFile);
  pathnames <- pathnames[keep];
  if (length(pathnames) == 0)
    throw("Failed to located ROC data set: ", tags);
  pathname <- pathnames[1];
  pathname;
}


loadLogRatios <- function(set, ...) {
  cat("Reading log-ratios: ", set$name, "\n", sep="");

  pathname <- Arguments$getReadablePathname(set$pathname);
  hasLogRatios <- (regexpr(",M,", pathname) != -1);

  if (hasLogRatios) {
    # Load log-ratios
    M <- loadObject(pathname);
  } else {
    # Load theta:s
    theta <- loadObject(pathname);

    # Get reference theta:s
    # We use females only, because that is what GTC does
    females <- (n23[colnames(theta)] == 2);
    thetaR <- matrixStats::rowMedians(theta[,females], na.rm=TRUE);

    # Calculate log-ratios and log-intensities
    M <- log2(theta/thetaR);

    rm(theta, thetaR);
  }
  rm(pathname);

  # Keep only subset to study
  M <- M[unitsToKeep,,drop=FALSE];

  # Match to samples of interest (in correct order)
  idxs <- match(sampleNames, colnames(M));
  if (anyMissing(idxs)) {
    print(sampleNames);
    print(colnames(M));
    throw("Internal error. Missing sample names: ", sampleNames[is.na(idxs)]);
  }
  M <- M[,idxs,drop=FALSE];
  rm(idxs);

  # Sanity checks
  if (ncol(M) != nbrOfSamples) {
    throw("Internal error. Incorrect number of samples in log-ratios: ", 
                                          ncol(M), " != " , nbrOfSamples);
  }
  if (nrow(M) != nbrOfUnits) {
    throw("Internal error. Incorrect number of units in log-ratios: ", 
                                             nrow(M), " != ", nbrOfUnits);
  }

  M;
} # loadLogRatios()
