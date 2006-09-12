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
#*/###########################################################################
setConstructorS3("SnpChipEffectFile", function(..., mergeStrands=FALSE) {
  extend(ChipEffectFile(...), "SnpChipEffectFile",
    mergeStrands = mergeStrands
  )
})


setMethodS3("getCellIndices", "SnpChipEffectFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If merging strands, return only every second group
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, function(groups) {
      groups[seq(from=1, to=length(groups), by=2)];
    })
  }

  cells;
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
