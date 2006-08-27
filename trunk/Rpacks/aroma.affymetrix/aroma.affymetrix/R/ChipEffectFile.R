###########################################################################/**
# @RdocClass ChipEffectFile
#
# @title "The ChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getChipEffects()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ChipEffectFile", function(..., model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(ParameterCelFile(...), "ChipEffectFile",
    "cached:.firstCells" = NULL,
    model = model
  )
})


setMethodS3("getFirstCellIndices", "ChipEffectFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose': 
 verbose <- Arguments$getVerbose(verbose);

  res <- this$.firstCells;
  if (is.null(res)) {
    stratifyBy <- switch(this$model, pm="pm");
    cdf <- getCdf(this);
    res <- getFirstCellIndices(cdf, units=NULL, ..., stratifyBy=stratifyBy, verbose=verbose);
    this$.firstCells <- res;
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, protected=TRUE)



setMethodS3("encodeUnitGroup", "ChipEffectFile", function(static, groupData, ...) {
  theta <- .subset2(groupData, "theta");
  ncells <- length(theta);
  stdvs <- rep(1, ncells);
  pixels <- rep(0, ncells);
  list(intensities=theta, stdvs=stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)




setMethodS3("decodeUnitGroup", "ChipEffectFile", function(static, groupData, ...) {
  res <- list();
  if (!is.null(groupData$intensities))
    res$theta <- groupData$intensities;
  if (!is.null(groupData$stdvs))
    res$stdvs <- groupData$stdvs;
  if (!is.null(groupData$pixels))
    res$pixels <- groupData$pixels;
  res;
}, static=TRUE, protected=TRUE)



setMethodS3("readUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- getFirstCellIndices(this, units=units, ...);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  stratifyBy <- switch(this$model, pm="pm");
  NextMethod("readUnits", this, cdf=cdf, ..., stratifyBy=stratifyBy);
});



setMethodS3("updateUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- getFirstCellIndices(this, units=units, ...);
  }

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
}, protected=TRUE);


############################################################################
# HISTORY:
# 2006-08-26
# o Created.  Have to store chip-effect estimates too.  Currently we use
#   the existing CEL/CDF structure for this, but those are unnecessarily
#   large for this.  Later this will be done in special CEL files with a
#   custom CDF file (possible virtual).  This requires that affxparser can
#   create empty CEL files from scratch, which is on the to-do list.
############################################################################
