###########################################################################/**
# @RdocClass ProbePreprocessing
#
# @title "The ProbePreprocessing class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a preprocessor that transforms
#  probe-level signals, typically intensities.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of @see "Preprocessing".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \details{
#   Subclasses must implement the \code{process()} method.
# }
#
# @author
#*/###########################################################################
setConstructorS3("ProbePreprocessing", function(...) {
  extend(Preprocessing(...), "ProbePreprocessing")
}, abstract=TRUE)


setMethodS3("getRootPath", "ProbePreprocessing", function(this, ...) {
  "probeData";
}, private=TRUE)



############################################################################
# HISTORY:
# 2006-12-08
# o Created.
############################################################################
