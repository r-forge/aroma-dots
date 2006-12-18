###########################################################################/**
# @RdocClass SnpChipEffectFile
#
# @title "The SnpChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectFile".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "SnpChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("SnpChipEffectFile", function(..., mergeStrands=FALSE) {
  extend(ChipEffectFile(...), "SnpChipEffectFile",
    "cached:.cellIndices" = NULL,
    mergeStrands = mergeStrands
  )
})


setMethodS3("as.character", "SnpChipEffectFile", function(this, ...) {
  s <- NextMethod("as.character", ...);
  s <- c(s, sprintf("Merge strands: %s", this$mergeStrands));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("getCellIndices", "SnpChipEffectFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- digest(list(...));
  res <- this$.cellIndices[[key]];
  if (!is.null(res)) {
    verbose && cat(verbose, "getCellIndices.SnpChipEffectFile(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get and restructure cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- NextMethod("getCellIndices", this, ..., verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge strands?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If merging strands, we only need half the number of chip-effect 
  # parameters per unit group.  Example:
  # a) mergeStrands=FALSE:
  #   Fit by strand and allele:        #groups=4, #chip effects=4
  #   (same but single-stranded SNP)   #groups=2, #chip effects=2
  # b) mergeStrands=TRUE:
  #   Merge strands, fit by allele:    #groups=4, #chip effects=2
  #   (same but single-stranded SNP)   #groups=2, #chip effects=2
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
#      groups[1:ceiling(ngroups/2)];
      groups[1:round((ngroups+1)/2)];
    })
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "getCellIndices.SnpChipEffectFile(): Updating cache");
  this$.cellIndices <- list();
  this$.cellIndices[[key]] <- cells;

  cells;
})


setMethodS3("readUnits", "SnpChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Check for cached data
  key <- digest(list(class=class(this), mergeStrands=this$mergeStrands, ...));
  res <- this$.readUnitsCache[[key]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.SnpChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", this, ..., force=TRUE, cache=FALSE, verbose=verbose);

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.SnpChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[key]] <- res;
  }

  res;
})


############################################################################
# HISTORY:
# 2006-12-18
# o BUG FIX: getCellIndices() would return a single group instead of two,
#   for a single-stranded SNP when mergeStrands=TRUE.  See for instance
#   unit SNP_A-1780520 in the Mapping250K_Nsp chip.
# 2006-11-28
# o Added readUnits() to override caching mechanism of superclasses.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated.  Now the names of the groups reflects the allele names as 
#   expected.
# 2006-09-11
# o Created.
############################################################################
