setMethodS3("extractThetaMatrix", "ChipEffectFile", function(this, units=NULL, groups=NULL, ..., verbose=FALSE) {
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
  # Import (theta1,theta2,...)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- extractDataFrame(this, units=units, verbose=less(verbose, 5));
  if (!is.null(groups)) {
    idxs <- which(data$group %in% groups);
    data <- data[idxs,,drop=FALSE];
    rm(idxs);
  } else {
    groups <- sort(unique(data$group));
  }
  verbose && cat(verbose, "Raw data:");
  verbose && str(verbose, data);

  nbrOfGroups <- length(groups);
  theta <- matrix(NA, nrow=nbrOfUnits, ncol=nbrOfGroups);
  for (gg in groups) {
    idxs <- which(data$group == gg);
    units <- data[idxs,"unit"];
    values <- data[idxs,ncol(data)];
    theta[units,gg] <- values;
  }
  rm(data);
  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, theta);

  theta;
})


setMethodS3("extractThetaMatrix", "SnpChipEffectFile", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  theta <- NextMethod("extractThetaMatrix", this, groups=groups, ...);

  theta;
})



setMethodS3("extractThetaMatrix", "CnChipEffectFile", function(this, groups=NULL, ...) {
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

  theta <- NextMethod("extractThetaMatrix", this, groups=groups, ...);

  theta;
})


############################################################################
# HISTORY:
# 2008-05-09
# o Created.
############################################################################
