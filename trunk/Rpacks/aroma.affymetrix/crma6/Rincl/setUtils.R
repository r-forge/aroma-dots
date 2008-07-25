# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
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

addRocDataSet <- function(sets, key, name=key, ...) {
  # Already added?
  if (name %in% names(sets))
    return(sets);

  cat("Adding data set: ", name, "\n", sep="");

  # Assert that it exists
  pathname <- getRocPathname(key);
  set <- list(name=name, pathname=pathname);
  sets[[name]] <- set;

  sets;
}

loadLogRatios <- function(set, ...) {
  cat("Getting log-ratios: ", set$name, "\n", sep="");

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


addRocData <- function(set, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  res <- set$rocData;
  if (!force && !is.null(res))
    return(set);

  verbose && enter(verbose, "Building RocData");
  verbose && cat(verbose, "Data set: ", set$name);

  # Get log-ratios
  M <- loadLogRatios(set);

  # Get true CNs
  C <- set$truth;
  if (is.null(C)) {
    verbose && enter(verbose, "Creating default truth (by columns)");
    C <- rep(1, times=ncol(M));
    C[(n23[colnames(M)] == 2)] <- 0;
    C <- as.integer(C);
    verbose && str(verbose, C);
    verbose && exit(verbose);
  }

  res <- RocData(truth=C, data=M);
  verbose && exit(verbose);

  set$rocData <- res;

  set;
} # addRocData()

extractRocData <- function(set, ...) {
  set <- addRocData(set, ...);
  set$rocData;
}


getSubsetDataSet <- function(set, subset=c("snp", "cn"), ..., verbose=FALSE) {
  # Argument 'subset':
  subset <- match.arg(subset);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting subset");
  verbose && cat(verbose, "Data set: ", set$name);
  verbose && cat(verbose, "Subset: ", subset);
  name <- sprintf("%s,subset=%s", set$name, subset);
  verbose && cat(verbose, "New data set: ", name);

  # Subset by unit type?
  types <- attr(unitTypes, "map");
  keep <- whichVector(unitTypes == types[subset]);

  # Get log-ratios
  set <- addRocData(set, ...);
  rocData <- extractSubset(set$rocData, rows=keep);
  set$rocData <- rocData;
  rm(rocData);

  set$name <- name;

  # Get positions
  if (!is.null(set$positions))
    set$positions <- set$positions[keep];
  rm(keep);

  verbose && exit(verbose);

  set;
} # getSubsetDataSet()


addSubsetDataSets <- function(sets, subset, ..., verbose=FALSE) {
  # Nothing to do?
  if (subset == "all")
    return(sets);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  rawSets <- sets[isRawSet(sets) & !isSubsetSet(sets)];
  for (kk in seq(along=rawSets)) {
    set <- rawSets[[kk]];
    key <- sprintf("%s,subset=%s", set$name, subset);
    if (key %in% names(sets))
      next;

    verbose && enter(verbose, sprintf("Subsetting data set %d", kk));
    set <- getSubsetDataSet(set, subset=subset, ..., verbose=verbose);

    sets[[key]] <- set;
    rm(set);

    verbose && exit(verbose);
  }

  sets;
} # addSubsetDataSets()





getSmoothDataSet <- function(set, h=2, ..., verbose=FALSE) {
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(2, Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Creating smoothed data set");
  verbose && cat(verbose, "Amount of smoothing: ", h);

  # Initial data ordered along the genome
  x <- set$positions;
  if (is.null(x))
    x <- positions;
  verbose && cat(verbose, "Positions: ");
  verbose && str(verbose, x);
  
  o <- order(x);
  x0 <- x[o];
  rm(x);
  verbose && cat(verbose, "Ordered positions: ");
  verbose && str(verbose, x0);

  set <- addRocData(set);
  rocData <- set$rocData;
  M <- getData(rocData, ordered=FALSE, complete=FALSE);
  M0 <- M[o,,drop=FALSE];
  verbose && cat(verbose, "Ordered log-ratios: ");
  verbose && str(verbose, M0);

  C0 <- NULL;
  if (!isTruthByColumns(rocData)) {
    C0 <- getTruth(rocData, ordered=FALSE, complete=FALSE);
    C0 <- C0[o,,drop=FALSE];
  }
  rm(o);

  idxs <- getBlockAverageMap(n=nrow(M0), h=h);
  hApprox <- attr(idxs, "hApprox");

##  verbose && enter(verbose, "Smoothing genome locations");
##  set$positions <- blockAvg(x0, idxs)
##  verbose && exit(verbose);
  rm(x0);

  set <- list(
    name = sprintf("%s,h=%.3f", set$name, h),
    h = h, hApprox = hApprox
  );

  if (!is.null(C0)) {
    verbose && enter(verbose, "Smoothing true CNs");
    set$truth <- blockAvg(C0, idxs);
    verbose && exit(verbose);
  }
  rm(C0);

  verbose && enter(verbose, "Smoothing log-ratios");
  set$M <- blockAvg(M0, idxs);
  rm(M0);
  verbose && exit(verbose);

  verbose && exit(verbose);

  set;
} # getSmoothDataSet()


addSmoothDataSets <- function(sets, hs=2:4, subset="all", ..., verbose=FALSE) {
  # Argument 'hs':
  hs <- Arguments$getDoubles(hs, range=c(1, Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Add raw sets for subset, if missing
  sets <- addSubsetDataSets(sets, subset=subset, verbose=verbose);

  rawSets <- sets[isRawSet(sets) & isSubsetSet(sets, subset=subset)];
  for (kk in seq(along=rawSets)) {
    verbose && enter(verbose, sprintf("Smoothing data set #%d", kk));
    set <- rawSets[[kk]];
    verbose && cat(verbose, "Name: ", set$name);

    for (h in hs) {
      hStr <- sprintf("h=%.3f", h);
      verbose && enter(verbose, hStr);
      key <- sprintf("%s,%s", set$name, hStr);
      if (key %in% names(sets)) {
        verbose && exit(verbose);
        next;
      }
      sets[[key]] <- getSmoothDataSet(set, h=h, ..., verbose=verbose);
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  sets;
} # addSmoothDataSets()





isSmoothedSet <- function(sets, ...) {
  (regexpr("h=", names(sets)) != -1);
}

isRawSet <- function(sets, ...) {
  (regexpr("h=", names(sets)) == -1);
}

isSubsetSet <- function(sets, subset=NULL, ...) {
  if (identical(subset, "all")) {
    pattern <- paste("subset=", sep="");
    (regexpr(pattern, names(sets)) == -1);
  } else {
    pattern <- paste("subset=", subset, sep="");
    (regexpr(pattern, names(sets)) != -1);
  }
}

isRawDataSet <- function(sets, ...) {
  (regexpr(",", names(sets)) == -1);
}

removeRawDataSets <- function(sets, ...) {
  sets[!isRawDataSet(sets)];
}


updateGraphics <- function(sets, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  hasPrefix <- function(name, prefix, ...) {
    (substring(name, 1, nchar(prefix)) == prefix);
  }

  hasAsterisk <- function(name, ...) {
    (regexpr("*", name, fixed=TRUE) != -1);
  }

  hasPlus <- function(name, ...) {
    (regexpr("+", name, fixed=TRUE) != -1);
  }

  for (kk in seq(along=sets)) {
    key <- names(sets)[kk];
    set <- sets[[kk]];
    name <- set$name;
    verbose && enter(verbose, sprintf("%d %s\n", kk, name));
  
    # Default colors
    col <- 0;
    lty <- 3;

    if (hasPrefix(name, "CRMA")) {
      col <- colors["CRMA"];
      if (hasPlus(name) & !hasAsterisk(name)) {
        lty <- 1;
      } else if (hasPlus(name) & hasAsterisk(name)) {
        lty <- 2;
      } else {
        lty <- 3;
      }
    } else if (hasPrefix(name, "dChip")) {
      col <- colors["dChip"];
      lty <- 1;
    } else if (hasPrefix(name, "APT")) {
      col <- colors["APT"];
      lty <- 1;
    } else if (hasPrefix(name, "GTC")) {
      col <- colors["GTC"];
      lty <- 1;
    }

    if (!identical(col, set$col)) {
      set$col <- col;
    }
    if (!identical(lty, set$lty)) {
      set$lty <- lty;
    }

    sets[[key]] <- set;

    rm(set, col, lty, name);
    verbose && exit(verbose);
  } # for (kk ...)

  sets;
} # updateGraphics()
