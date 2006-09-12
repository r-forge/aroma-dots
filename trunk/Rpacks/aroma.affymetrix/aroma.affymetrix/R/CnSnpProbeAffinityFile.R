###########################################################################/**
# @RdocClass CnSnpProbeAffinityFile
#
# @title "The CnSnpProbeAffinityFile class"
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
setConstructorS3("SnpProbeAffinityFile", function(..., averageAB=FALSE) {
  extend(SnpProbeAffinityFile(...), "CnSnpProbeAffinityFile",
    averageAB=averageAB
  )
})


setMethodS3("getCellIndices", "SnpProbeAffinityFile", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If merging strands, return only every second group
  if (this$averageAB) {
    cells <- applyCdfGroups(cells, function(groups) {
      groups[seq(from=1, to=length(groups), by=2)];
    })
  }

  cells;
})

setMethodS3("setAverageAB", "SnpProbeAffinityFile", function(this, status, ...) {
  this$averageAB <- status;
})



############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
