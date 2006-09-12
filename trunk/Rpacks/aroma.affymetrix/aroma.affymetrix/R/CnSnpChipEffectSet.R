###########################################################################/**
# @RdocClass CnSnpChipEffectSet
#
# @title "The CnSnpChipEffectSet class"
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
#   \item{averageAB}{A @logical indicating if the signals from allele A and
#      allele B are averaged or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnSnpChipEffectSet", function(..., averageAB=FALSE) {
  this <- extend(SnpChipEffectSet(...), "CnSnpChipEffectSet");
  setAverageAB(this, averageAB);
  this;
})

setMethodS3("getChipEffectFileClass", "CnSnpChipEffectSet", function(static, ...) {
  CnSnpChipEffectFile;
}, static=TRUE)


setMethodS3("getAverageAB", "CnSnpChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$averageAB;
})

setMethodS3("setAverageAB", "CnSnpChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  status <- Arguments$getLogical(status);
  oldStatus <- getAverageAB(this);
  ce <- getFile(this, 1);
  ce$averageAB <- status;
  invisible(oldStatus);
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
