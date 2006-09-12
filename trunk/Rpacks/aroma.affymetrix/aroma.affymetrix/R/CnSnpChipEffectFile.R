###########################################################################/**
# @RdocClass CnSnpChipEffectFile
#
# @title "The CnSnpChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in a copy-number probe-level
#  models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpChipEffectFile".}
#   \item{averageAB}{A @logical indicating if the signals from allele A and
#      allele B are averaged or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnSnpChipEffectFile", function(..., averageAB=FALSE) {
  extend(SnpChipEffectFile(...), "CnSnpChipEffectFile",
    averageAB = averageAB
  )
})

setMethodS3("getCellIndices", "CnSnpChipEffectFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If merging strands, return only every second group
  if (this$averageAB) {
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
