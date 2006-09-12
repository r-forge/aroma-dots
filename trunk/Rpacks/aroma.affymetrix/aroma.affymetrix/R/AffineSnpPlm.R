###########################################################################/**
# @RdocClass AffineSnpPlm
#
# @title "The AffineSnpPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffinePlm".}
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
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "AffineProbeAffinityFile" class.  
# }
#
# \references{
# }
#*/###########################################################################
setConstructorS3("AffineSnpPlm", function(..., name="modelAffineSnpPlm", mergeStrands=FALSE) {
  extend(AffinePlm(..., name=name), c("AffineSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created from the MbeiSnpPlm.
############################################################################
