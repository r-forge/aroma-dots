setMethodS3("extractTheta", "ChipEffectFile", function(this, units=NULL, groups=NULL, ..., verbose=FALSE) {
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
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'groups':
  if (!is.null(groups)) {
    groups <- Arguments$getIndices(groups, range=c(1, 999));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }


  # Subset by groups?
  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Import (theta1,theta2,...)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- extractMatrix(this, units=ugcMap, verbose=less(verbose, 5));
  data <- data[,1,drop=TRUE];
  verbose && cat(verbose, "Raw data:");
  verbose && str(verbose, data);

  allUnits <- unique(ugcMap[,"unit"]);
  nbrOfGroups <- length(groups);
  theta <- matrix(NA, nrow=nbrOfUnits, ncol=nbrOfGroups);
  for (gg in groups) {
    idxs <- which(ugcMap$group == gg);
    units <- ugcMap[idxs,"unit"];
    units <- match(units, allUnits);
    theta[units,gg] <- data[idxs];
    rm(idxs, units);
  }
  rm(data, allUnits);

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, theta);

  theta;
})




setMethodS3("extractTheta", "SnpChipEffectFile", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  theta <- NextMethod("extractTheta", this, groups=groups, ...);

  theta;
})



setMethodS3("extractTheta", "CnChipEffectFile", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    if (this$combineAlleles) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  theta <- NextMethod("extractTheta", this, groups=groups, ...);

  theta;
})



############################################################################
# HISTORY:
# 2008-05-10
# o Updated to take an UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
