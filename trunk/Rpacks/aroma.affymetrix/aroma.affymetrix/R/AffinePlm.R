###########################################################################/**
# @RdocClass AffinePlm
#
# @title "The AffinePlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Bengtsson \& Hössjer (2006) model,
#  see @see "AffinePlm".
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
setConstructorS3("AffinePlm", function(..., name="modelAffine") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(aroma.light) || throw("Package 'aroma.light' not loaded.");

  extend(ProbeLevelModel(..., name=name), "AffinePlm")
})


setMethodS3("getProbeAffinityClass", "AffinePlm", function(static, ...) {
  AffineProbeAffinityFile;
}, static=TRUE, protected=TRUE)


setMethodS3("getFitFunction", "AffinePlm", function(static, ...) {
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
