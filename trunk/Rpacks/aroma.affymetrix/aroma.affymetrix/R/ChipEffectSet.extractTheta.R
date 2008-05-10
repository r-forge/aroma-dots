setMethodS3("extractTheta", "ChipEffectSet", function(this, units=NULL, groups=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);
  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
    nbrOfUnits <- length(units);
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
  ugcMap <- NULL;
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

  verbose && cat(verbose, "Filtered (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);


  nbrOfGroups <- length(groups);
  nbrOfArrays <- nbrOfArrays(this);
  dim <- c(nbrOfUnits, nbrOfGroups, nbrOfArrays);
  dimnames <- list(NULL, NULL, getNames(this));
  theta <- array(NA, dim=dim, dimnames=dimnames);
  for (kk in seq(length=nbrOfArrays)) {
    ce <- getFile(this, kk);
    thetaKK <- extractTheta(ce, units=ugcMap, groups=groups, verbose=less(verbose, 5));
    verbose && str(verbose, thetaKK);
    theta[,,kk] <- thetaKK;
  }
  rm(ugcMap);

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, theta);

  theta;
})



setMethodS3("extractTheta", "SnpChipEffectSet", function(this, groups=NULL, ...) {
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



setMethodS3("extractTheta", "CnChipEffectSet", function(this, groups=NULL, ...) {
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
# o Created.
############################################################################
