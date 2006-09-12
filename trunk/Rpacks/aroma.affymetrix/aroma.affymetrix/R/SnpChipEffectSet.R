###########################################################################/**
# @RdocClass SnpChipEffectSet
#
# @title "The SnpChipEffectSet class"
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
#   \item{...}{Arguments passed to @see "ChipEffectSet".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("SnpChipEffectSet", function(..., mergeStrands=FALSE) {
  this <- extend(ChipEffectSet(...), "SnpChipEffectSet");
  setMergeStrands(this, mergeStrands);
  this;
})


setMethodS3("getChipEffectFileClass", "SnpChipEffectSet", function(static, ...) {
  SnpChipEffectFile;
}, static=TRUE)

setMethodS3("getMergeStrands", "SnpChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$mergeStrands;
})

setMethodS3("setMergeStrands", "SnpChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);

  # Argument 'status':
  status <- Arguments$getLogical(status);

  oldStatus <- getMergeStrands(this);

  # Update all chip-effect files
  lapply(this, function(ce) {
    ce$mergeStrands <- status;
  })

  invisible(oldStatus);
})



############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
