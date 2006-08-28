###########################################################################/**
# @RdocClass AffymetrixAffineModel
#
# @title "The AffymetrixAffineModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Bengtsson \& Hössjer (2006) model,
#  see @see "AffymetrixAffineModel".
# }
# 
# @synopsis
#
# \arguments{
#   \item{path}{The @character string specifying the path to the directory
#      to contain the parameter-estimate files.}
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \references{
#   Bengtsson \& Hössjer (2006). \cr
# }
#*/###########################################################################
setConstructorS3("AffymetrixAffineModel", function(..., name="modelAffine") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(aroma.light) || throw("Package 'aroma.light' not loaded.");

  extend(ProbeLevelModel(..., name=name), "AffymetrixAffineModel")
})


setMethodS3("getProbeAffinityClass", "AffymetrixAffineModel", function(static, ...) {
  AffineProbeAffinityFile;
}, static=TRUE, protected=TRUE)


setMethodS3("getFitFunction", "AffymetrixAffineModel", function(static, ...) {
  fcn <- function(y, ...) {
    f <- calibrateMultiscan(t(y), center=FALSE, project=TRUE);
    phi <- as.vector(attr(f, "modelFit")$b);
    list(
      theta = as.vector(f), 
      phi = phi
    )
  }

  fcn;
}, static=TRUE, protected=TRUE)





############################################################################
# HISTORY:
# 2006-08-28
# o Created from the Li & Wong model.
############################################################################
