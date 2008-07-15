setMethodS3("getUnitGroupCellMatrixMap", "ChipEffectFile", function(this, units=NULL, groups=NULL, ..., verbose=FALSE) {
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
  # Build integer UxG matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  allUnits <- unique(ugcMap[,"unit"]);
  nbrOfGroups <- length(groups);
  map <- matrix(as.integer(NA), nrow=nbrOfUnits, ncol=nbrOfGroups);
  
  for (gg in seq(length=nbrOfGroups)) {
    group <- groups[gg];
    verbose && enter(verbose, sprintf("Group %d (%d) of %d", 
                                                  gg, group, nbrOfGroups));

    idxs <- which(ugcMap$group == group);
    units <- ugcMap[idxs, "unit"];
    cells <- ugcMap[idxs, "cell"];
    rr <- match(units, allUnits);
    map[rr,gg] <- cells;

    rm(idxs, rr, units, cells);
    verbose && exit(verbose);
  }

  class(map) <- "UnitGroupCellMatrixMap";

  verbose && cat(verbose, "Unit-by-group cell matrix map:");
  verbose && str(verbose, map);

  map;
}, protected=TRUE)  # getUnitGroupCellMatrixMap()




setMethodS3("extractTheta", "ChipEffectFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'units':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (inherits(units, "UnitGroupCellMatrixMap")) {
    cellMatrixMap <- units;
  } else {
    cellMatrixMap <- getUnitGroupCellMatrixMap(this, units=units, ..., 
                                               verbose=less(verbose, 10));
  }

  data <- readRawData(this, indices=cellMatrixMap, fields="intensities", 
                                    drop=TRUE, verbose=less(verbose, 20));
  dim(data) <- dim(cellMatrixMap);

  # Drop singleton dimensions
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, data);

  data;
}) # extractTheta()




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
}) # extractTheta()



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
}) # extractTheta()



############################################################################
# HISTORY:
# 2008-07-13
# o Added argument 'drop=FALSE' to extractTheta().
# 2008-06-09
# o Added getUnitGroupCellMatrixMap() to ChipEffectFile.  The extractTheta()
#   methods is now using this method.
# 2008-05-10
# o Updated to take an UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
