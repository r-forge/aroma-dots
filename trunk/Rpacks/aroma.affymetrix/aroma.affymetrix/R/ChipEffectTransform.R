###########################################################################/**
# @RdocClass ChipEffectTransform
#
# @title "The ChipEffectTransform class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a transform that transforms chip-effect
#  estimates obtained from probe-level modelling.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{The input data set as an @see "ChipEffectSet".}
#   \item{...}{Arguments passed to the constructor of @see "Transform".}
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
setConstructorS3("ChipEffectTransform", function(dataSet=NULL, ...) {
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

  extend(Transform(dataSet=dataSet, ...), "ChipEffectTransform")
}, abstract=TRUE)


setMethodS3("getRootPath", "ChipEffectTransform", function(this, ...) {
  "plmData";
}, private=TRUE)


setMethodS3("getOutputDataSet", "ChipEffectTransform", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting output data set");

  args <- list(generic="getOutputDataSet", this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inherit certain arguments from the input data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AD HOC (not using OO), but setting these arguments does speed
  # up things. /HB 2007-09-17
  # Note, this is also done in Transform for now, such it is not really
  # needed here.  However, in case it will be removed from there it still
  # makes sense to have it here.
  ds <- getInputDataSet(this);
  if (inherits(ds, "CnChipEffectSet"))
    args$combineAlleles <- ds$combineAlleles;
  if (inherits(ds, "SnpChipEffectSet"))
    args$mergeStrands <- ds$mergeStrands; 

  verbose && cat(verbose, "Calling NextMethod() with arguments:");
  verbose && str(verbose, args);

  args$verbose <- less(verbose, 10);
  res <- do.call("NextMethod", args);

  verbose && exit(verbose);

  outputDataSet;
})


############################################################################
# HISTORY:
# 2007-09-18
# o Now getOutputDataSet() of Transform carry down certain arguments from
#   the input data set. This will speed up things.
# 2006-12-08
# o Created.
############################################################################
