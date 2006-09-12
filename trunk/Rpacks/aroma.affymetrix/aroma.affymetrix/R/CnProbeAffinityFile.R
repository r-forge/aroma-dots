###########################################################################/**
# @RdocClass CnProbeAffinityFile
#
# @title "The CnProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in SNP probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpProbeAffinityFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnProbeAffinityFile", function(..., combineAlleles=FALSE) {
  extend(SnpProbeAffinityFile(...), "CnProbeAffinityFile",
    combineAlleles=combineAlleles
  )
})


setMethodS3("getCellIndices", "CnProbeAffinityFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If merging strands, return only every second group
  if (this$combineAlleles) {
    cells <- applyCdfGroups(cells, function(groups) {
      groups[seq(from=1, to=length(groups), by=2)];
    })
  }

  cells;
})

setMethodS3("setCombineAlleles", "CnProbeAffinityFile", function(this, status, ...) {
  this$combineAlleles <- status;
})



############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
