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

  # If combining alleles, return only every second group.
  # In order to improve readability we merge the names of alleles groups
  # combined, e.g. groups 'C' and 'G' become group 'CG'.
  if (this$combineAlleles) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      odds <- seq(from=1, to=ngroups, by=2);
      evens <- seq(from=2, to=ngroups, by=2);
      names <- names(groups);
      names <- paste(names[odds], names[evens], sep="");
      groups <- groups[odds];
      names(groups) <- names;
      groups;
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
