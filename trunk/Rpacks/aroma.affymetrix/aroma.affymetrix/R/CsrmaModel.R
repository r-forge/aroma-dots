###########################################################################/**
# @RdocClass CsrmaModel
#
# @title "The CsrmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Chromosomal Smoothing Robust Multichip Analysis
#  method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#*/###########################################################################
setConstructorS3("CsrmaModel", function(cesTuple=NULL, ...) {
  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "CsrmaModel")
})


setMethodS3("getAsteriskTag", "CsrmaModel", function(this, ...) {
  "CSRMA";
}, protected=TRUE)



##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
