###########################################################################/**
# @RdocClass SmoothRmaModel
#
# @title "The SmoothRmaModel class"
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
#   \item{...}{Arguments passed to the constructor of 
#              @see "SmoothMultiarrayModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("SmoothRmaModel", function(...) {
  extend(SmoothMultiarrayModel(...), "SmoothRmaModel");
})


setMethodS3("getAsteriskTag", "SmoothRmaModel", function(this, ...) {
  tags <- NextMethod("getAsteriskTag", this, ...);
  tags[1] <- "SRMA";
  tags;
}, protected=TRUE)


setMethodS3("getFitFunction", "SmoothRmaModel", function(this, ...) {
  smoothWRMA;
})


##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
