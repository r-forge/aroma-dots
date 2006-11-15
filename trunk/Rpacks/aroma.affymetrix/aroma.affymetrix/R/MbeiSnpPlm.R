###########################################################################/**
# @RdocClass MbeiSnpPlm
#
# @title "The MbeiSnpPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiPlm".}
#   \item{name}{The name of the model, which is also used in the pathname.}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("MbeiSnpPlm", function(..., name="modelMbeiSnpPlm", mergeStrands=FALSE) {
  extend(MbeiPlm(..., name=name), c("MbeiSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})

setMethodS3("getRootPath", "MbeiSnpPlm", function(this, ...) {
  "modelMbeiSnpPlm";
})

############################################################################
# HISTORY:
# 2006-09-11
# o Created from the RmaSnpPlm.
############################################################################
