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

  # If merging strands, we only need half the number of chip-effect 
  # parameters per unit group.
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      groups[1:ceiling(ngroups/2)];
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


############################################################################
# HISTORY:
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated.  Now the names of the groups reflects the allele names as 
#   expected.
# 2006-09-11
# o Created.
############################################################################
