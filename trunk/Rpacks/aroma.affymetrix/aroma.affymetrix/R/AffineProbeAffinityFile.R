###########################################################################/**
# @RdocClass AffineProbeAffinityFile
#
# @title "The AffineProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in the 
#  Bengtsson \& Hössjer (2006) model.
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
#   An object of this class is typically obtain through the
#   \code{getProbeAffinities()} method for the 
#   @see "AffymetrixAffineModel" class.
# }
#
#*/###########################################################################
setConstructorS3("AffineProbeAffinityFile", function(...) {
  extend(ProbeAffinityFile(...), "AffineProbeAffinityFile");
})


setMethodS3("encodeUnitGroup", "AffineProbeAffinityFile", function(static, groupData, ...) {
  phi <- .subset2(groupData, "phi");
  ncells <- length(phi);
  stdvs <- rep(1, ncells);
  pixels <- rep(0, ncells);
  list(intensities=phi, stdvs=stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)




setMethodS3("decodeUnitGroup", "AffineProbeAffinityFile", function(static, groupData, ...) {
  res <- list();
  if (!is.null(groupData$intensities))
    res$phi <- groupData$intensities;
  if (!is.null(groupData$stdvs))
    res$stdvs <- groupData$stdvs;
  if (!is.null(groupData$pixels))
    res$pixels <- groupData$pixels;
  # Rescale
  phi <- res$phi;
  theta <- res$theta;
  scale <- prod(phi);
  phi <- phi/scale;
  theta <- theta * scale;
  res$phi <- phi;
  res$theta <- theta;
  res;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2006-08-28
# o Created from RmaProbeAffinityFile.
############################################################################
