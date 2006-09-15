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
#   An object of this class is typically part of a @see "SnpChipEffectsSet".
# }
#*/###########################################################################
setConstructorS3("SnpChipEffectFile", function(..., mergeStrands=FALSE) {
  extend(ChipEffectFile(...), "SnpChipEffectFile",
    mergeStrands = mergeStrands
  )
})


setMethodS3("getCellIndices", "SnpChipEffectFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If merging strands, we only need half the number of chip-effect 
  # parameters per unit group.
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      groups[1:ceiling(ngroups/2)];
    })
  }

  cells;
})


############################################################################
# HISTORY:
# 2006-09-12
# o Updated.  Now the names of the groups reflects the allele names as 
#   expected.
# 2006-09-11
# o Created.
############################################################################
