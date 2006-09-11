###########################################################################/**
# @RdocClass MbeiPlm
#
# @title "The MbeiPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Li \& Wong (2001) model,
#  see @see "MbeiPlm".
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
#   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
#   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
# }
#*/###########################################################################
setConstructorS3("MbeiPlm", function(..., name="modelMbei") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affy) || throw("Package 'affy' not loaded.");

  extend(ProbeLevelModel(..., name=name), "MbeiPlm")
})


###########################################################################/**
# @RdocMethod getProbeAffinityClass
#
# @title "Static method to get the ProbeAffinityFile Class object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "Class" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeAffinityClass", "MbeiPlm", function(static, ...) {
  MbeiProbeAffinityFile;
}, static=TRUE, protected=TRUE)



###########################################################################/**
# @RdocMethod getFitFunction
#
# @title "Static method to get the low-level function that fits the PLM"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @function.
# }
#
# @author
#
# \seealso{
#   @see "affy::fit.li.wong".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "MbeiPlm", function(static, ...) {
  liWong <- function(y, ...) {
    fit <- fit.li.wong(t(y));

    # A fit function must return: theta, sdTheta, thetaOutliers, phi, sdPhi, phiOutliers.
    names <- names(fit);
    idxs <- match(c("sigma.theta", "theta.outliers", "sigma.phi", "phi.outliers"), names);
    names[idxs] <- c("sdTheta", "thetaOutliers", "sdPhi", "phiOutliers");
    names(fit) <- names;

    fit;
  }

  liWong;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2006-09-10
# o Updated getFitFunction() to return required fields.
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added getProbeAffinities() and the corrsponding cached fields.
# o Now fit() does not re-read data just updated.
# 2006-08-19
# o After all the bug fixes in updateCel() I think this function finally
#   works correctly.
# 2006-08-17
# o Created.
############################################################################
