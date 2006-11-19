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
setConstructorS3("MbeiSnpPlm", function(..., tags=c("MBEI", ifelse(mergeStrands, "", "+-")), mergeStrands=FALSE) {
  extend(MbeiPlm(..., tags=tags), c("MbeiSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created from the RmaSnpPlm.
############################################################################
