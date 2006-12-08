###########################################################################/**
# @RdocClass ChipEffectPreprocessing
#
# @title "The ChipEffectPreprocessing class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a preprocessor that transforms
#  chip effects obtained from probe-level modelling.
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
setConstructorS3("ChipEffectPreprocessing", function(dataSet=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "ChipEffectSet")) {
      throw("Argument 'dataSet' is not a ChipEffectSet object: ", 
                                                           class(dataSet));
    }
  }

  extend(Preprocessing(dataSet=dataSet, ...), "ChipEffectPreprocessing")
}, abstract=TRUE)


setMethodS3("getRootPath", "ChipEffectPreprocessing", function(this, ...) {
  "plmData";
})



############################################################################
# HISTORY:
# 2006-12-08
# o Created.
############################################################################
