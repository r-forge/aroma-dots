###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod isResequenceChip
#
# @title "Static method to check if a chip is a resequence (CustomSeq) chip"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the chip type refers to a resequence array, 
#   otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isResequenceChip", "AffymetrixCdfFile", function(this, ...) {
  chipType <- getChipType(this);

  # First some hardwired return values
  if (regexpr("^Mitochip_2.*$", chipType) != -1)
    return(TRUE);

  # Then, check for resequencing units
  types <- getUnitTypes(this, ...);
  hasReseqUnits <- any(types == 3);

  hasReseqUnits;
}, private=TRUE)





###########################################################################/**
# @RdocMethod readUnitsByQuartets
#
# @title "Gets the cell quartets for each base position"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{Subset of units to be queried. If @NULL, all units are used.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @vector of @factors.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readUnitsByQuartets", "AffymetrixCdfFile", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(this)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Read (pbase, tbase, cell index)
  verbose && enter(verbose, "Reading (pbase, tbase, cell index)");
  cdfUnits <- readUnits(this, units=units, readIndices=TRUE, readBases=TRUE, readXY=FALSE, readExpos=FALSE, readType=FALSE, readDirection=FALSE);
  verbose && cat(verbose, "Number of units read: ", length(cdfUnits));
  verbose && exit(verbose);

  verbose && enter(verbose, "Restructuring (and validating assumptions about) fields 'pbase', 'tbase', and 'indices'");
  verbose && cat(verbose, "Number of units: ", length(cdfUnits));
  # Restructure and validate
  for (uu in seq(along=cdfUnits)) {
    verbose && enter(verbose, sprintf("Unit #%d ('%s') of %d", uu, names(cdfUnits)[uu], length(cdfUnits)));

    cdfUnit <- cdfUnits[[uu]];
    cdfGroups <- cdfUnit$groups;
    verbose && cat(verbose, "Number of groups: ", length(cdfGroups));

    for (gg in seq(along=cdfGroups)) {
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d", gg, names(cdfGroups)[gg], length(cdfGroups)));

      cdfGroup <- cdfGroups[[gg]];
      # Sanity check of assumption of ordering of cells
      pbase <- cdfGroup$pbase;
      pbase <- matrix(pbase, nrow=4, byrow=FALSE);
      pbase <- t(pbase);
      pbase <- unique(pbase);
      if (nrow(pbase) != 1) {
        throw("Assumption exception: The probes are not ordered consistently.");
      }
      pbase <- as.vector(pbase);

      tbase <- cdfGroup$tbase;
      tbase <- matrix(tbase, nrow=4, byrow=FALSE);
      tbase <- tbase[1,,drop=TRUE];

      cells <- cdfGroup$indices;
      cells <- matrix(cells, nrow=4, byrow=FALSE);
      rownames(cells) <- pbase;

      cdfGroup <- list(indices=cells, tbase=tbase);

      cdfGroups[[gg]] <- cdfGroup;
      rm(pbase, cells, cdfGroup);
      verbose && exit(verbose);
    } # for (gg ...)

    cdfUnit$groups <- cdfGroups;
    cdfUnits[[uu]] <- cdfUnit;
    rm(cdfUnit);

    verbose && exit(verbose);
  } # for (uu ...)
  verbose && exit(verbose);

  cdfUnits;
}, private=TRUE)



setMethodS3("getCellQuartets", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(this)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting cell-index matrix");

  key <- list(method="getCellQuartets", class=class(this)[1], units=units);
  dirs <- c("aroma.affymetrix", getChipType(this));
  if (!force) {
    cells <- loadCache(key=key, dirs=dirs);
    if (!is.null(cells)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(cells);
    }
  }

  cdfUnits <- readUnitsByQuartets(this, units=units, verbose=verbose);

  pbase <- rownames(cdfUnits[[1]]$groups[[1]]$indices);
  # Sanity check
  if (is.null(pbase)) {
    throw("No resequencing cell indices available.");
  }

  verbose && enter(verbose, "Restructuring into a matrix");

  tbase <- applyCdfGroups(cdfUnits, cdfGetFields, "tbase");
  tbase <- unlist(tbase, use.names=FALSE);
  nbrOfBases <- length(tbase);

  cells <- applyCdfGroups(cdfUnits, cdfGetFields, "indices");
  rm(cdfUnits);  # Not needed anymore

  cells <- unlist(cells, use.names=FALSE);
  dim(cells) <- c(length(pbase), nbrOfBases);
  dimnames(cells) <- list(pbase, tbase);

  # Transpose
  cells <- t(cells);

  verbose && str(verbose, cells);
  verbose && exit(verbose);


  # Save to cache?
  if (cache) {
    saveCache(cells, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  cells;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-08-18
# BUG FIX: Used non-existing 'cdf' instead of 'this'.
# 2008-08-10
# o Created.
############################################################################
