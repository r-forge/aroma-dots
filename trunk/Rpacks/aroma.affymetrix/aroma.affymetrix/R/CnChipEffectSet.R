###########################################################################/**
# @RdocClass CnChipEffectSet
#
# @title "The CnChipEffectSet class"
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
#   \item{...}{Arguments passed to @see "SnpChipEffectSet".}
#   \item{combineAlleles}{A @logical indicating if the signals from allele A and
#      allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnChipEffectSet", function(..., combineAlleles=FALSE) {
  this <- extend(SnpChipEffectSet(...), "CnChipEffectSet");
  setCombineAlleles(this, combineAlleles);
  this;
})

setMethodS3("getChipEffectFileClass", "CnChipEffectSet", function(static, ...) {
  CnChipEffectFile;
}, static=TRUE)


setMethodS3("getCombineAlleles", "CnChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$combineAlleles;
})

setMethodS3("setCombineAlleles", "CnChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  status <- Arguments$getLogical(status);
  ce <- getFile(this, 1);
  oldStatus <- ce$ombineAlleles;
  ce$combineAlleles <- status;
  invisible(oldStatus);
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
