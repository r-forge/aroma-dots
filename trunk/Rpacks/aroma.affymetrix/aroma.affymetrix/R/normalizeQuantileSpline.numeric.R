###########################################################################/**
# @set "class=numeric"
# @RdocMethod normalizeQuantileSpline
#
# @title "Normalizes the empirical distribution of a single sample to a target distribution"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{x}{a @numeric @vector of length \eqn{N}.}
#   \item{xTarget}{a @numeric @vector of length \eqn{M}.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length \eqn{N}.
# }
#
# \section{Missing values}{
#   Both argument \code{X} and \code{xTarget} may contain non-finite values.
#   These values do not affect the estimation of the normalization function.
#   Non-finite values in \code{X}, remain in the output.
# }
#
# \seealso{
#   @see "normalizeQuantileSpline.matrix".
# }
#
# @author
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("normalizeQuantileSpline", "numeric", function(x, xTarget, sort=TRUE, ..., robust=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- length(x);

  # Argument 'xTarget':
  if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }

  if (length(xTarget) != n) {
    throw("Argument 'xTarget' is of different length than 'x': ", 
                                               length(xTarget), " != ", n);
  }

  # Sort signals to be normalized
  xx <- sort(x, na.last=TRUE);

  # Sort target distribution?
  if (sort)
    xTarget <- sort(xTarget, na.last=TRUE);

  # Keep only finite values
  ok <- (is.finite(xx) & is.finite(xTarget));
  xx <- xx[ok];
  xTarget <- xTarget[ok];
  rm(ok); # Not needed anymore

  if (robust) {
    fit <- smooth.spline(x=xx, y=xTarget, ...);
  } else {
    fit <- robustSmoothSpline(x=xx, y=xTarget, ...);
  }

  # Not needed anymore
  rm(xx, xTarget);

  # Normalize the data
  ok <- is.finite(x);
  x[ok] <- predict(fit, x=x[ok])$y;

  x;
}) # normalizeQuantileSpline.numeric()



##############################################################################
# HISTORY:
# 2007-02-05
# o Now normalizeQuantileSpline() handles NAs too.
# 2007-02-04
# o Created from normalizeQuantile.numeric.R.
##############################################################################
