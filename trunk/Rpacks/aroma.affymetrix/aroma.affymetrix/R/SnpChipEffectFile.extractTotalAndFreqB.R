setMethodS3("extractTotalAndFreqB", "CnChipEffectFile", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    rm(units);
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (this$combineAlleles && this$mergeStrands) {
    # theta == (theta)
    groups <- 1;
  } else if (this$combineAlleles && !this$mergeStrands) {
    # theta == (theta+, theta-)
    groups <- c(1,3);
  } else if (!this$combineAlleles && this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else if (!this$combineAlleles && !this$mergeStrands) {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  theta <- extractTheta(this, units=ugcMap, groups=groups, ..., 
                                                 verbose=less(verbose, 5));
  nbrOfUnits <- nrow(theta);

  # Calculating total chip effect
  thetaTotal <- rowSums(theta, na.rm=TRUE);

  # Calculating Allele B frequency
  if (this$combineAlleles) {
    freqB <- rep(NA, nbrOfUnits);
    rm(theta);
  } else {
    if (ncol(theta) == 2) {
      thetaB <- theta[,2];
    } else if (ncol(theta) == 4) {
      thetaB <- rowSums(theta[,c(2,4)], na.rm=TRUE);
    }
    rm(theta);
    freqB <- thetaB/thetaTotal;
    rm(thetaB);
  }

  data <- matrix(c(thetaTotal, freqB), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("total", "freqB");

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})




setMethodS3("extractTotalAndFreqB", "SnpChipEffectFile", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    rm(units);
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  theta <- extractTheta(this, units=ugcMap, groups=groups, ..., 
                                                 verbose=less(verbose, 5));
  nbrOfUnits <- nrow(theta);

  # Calculating total chip effect
  thetaTotal <- rowSums(theta, na.rm=TRUE);

  # Calculating Allele B frequencies
  if (ncol(theta) == 2) {
    thetaB <- theta[,2];
  } else if (ncol(theta) == 4) {
    thetaB <- rowSums(theta[,c(2,4)], na.rm=TRUE);
  }
  rm(theta);
  freqB <- thetaB/thetaTotal;
  rm(thetaB);

  data <- matrix(c(thetaTotal, freqB), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("total", "freqB");

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})





############################################################################
# HISTORY:
# 2008-05-10
# o Now extractTotalAndFreqB() takes and UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
