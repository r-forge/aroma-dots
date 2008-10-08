#########################################################################/**
# @set "class=matrix"
# @RdocMethod fitPrincipalCurve
#
# \encoding{latin1}
#
# @title "Fit a principal curve in K dimensions"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix (K>=2) where the columns represent the dimension.}
#  \item{fixDimension}{An optional @integer specifying which dimension
#    to keep fix such that the corresponding component of the estimated 
#    principal curve is (approximately) the identify function.}
#  \item{...}{Other arguments passed to @see "princurve::principal.curve".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a principal.curve object (which is a @list).
#   See @see "princurve::principal.curve" for more details.
# }
#
# \section{Missing values}{
#  The estimation of the affine normalization function will only be made
#  based on complete observations, i.e. observations that contains no @NA
#  values in any of the channels.
# }
#
# @author
#
# \references{
#   [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
# }
#
# @examples "../incl/fitPrincipalCurve.matrix.Rex"
#
# \seealso{
#   @seemethod "backtransformPrincipalCurve".
#   @see "princurve::principal.curve".
# }
#*/######################################################################### 
setMethodS3("fitPrincipalCurve", "matrix", function(X, fixDimension=NULL, ..., verbose=FALSE) {
  require("princurve") || throw("Package not loaded: princurve");

  # The current implementation contains bugs. /HB 2008-05-26
  principal.curve <- principal.curve.hb;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  n <- nrow(X);
  p <- ncol(X);

  if (!is.null(fixDimension)) {
    fixDimension <- Arguments$getIndex(fixDimension, range=c(1,p));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting principal curve");
  verbose && cat(verbose, "Data size: ", n, "x", p);

  if (!is.null(fixDimension)) {
    verbose && cat(verbose, "Keeping dimension fix: ", fixDimension);
  }

  verbose && enter(verbose, "Identifying missing values");
  # princurve::principal.curve() does not handle missing values.
  keep <- rep(TRUE, n);
  for (cc in seq(length=p))
    keep <- keep & is.finite(X[,cc]);
  anyMissing <- (!all(keep));
  if (anyMissing)
    X <- X[keep,];
  verbose && exit(verbose);

  verbose && cat(verbose, "Data size after removing non-finite data points: ", nrow(X), "x", p);


  verbose && enter(verbose, "Calling principal.curve()");
  fit <- principal.curve(X, ...);
  verbose && exit(verbose);

  fit$fixDimension <- fixDimension;

  if (!is.null(fixDimension)) {
    verbose && enter(verbose, "Fixing one of the dimensions");
    verbose && cat(verbose, "Dimension: ", fixDimension);

    verbose && enter(verbose, "Fitting 'fix' dimension");
    lambda <- fit$lambda;
    x <- fit$s[,fixDimension, drop=TRUE];
    if (length(x) != length(lambda)) {
      throw("Error in assumption. Internal 'x' and fitted 'lambda' are of different lengths: ", length(x), " != ", length(lambda));
    }
    ok <- whichVector(is.finite(x) & is.finite(lambda));
    lambda <- lambda[ok];
    x <- x[ok];
    rm(ok);
    if (length(x) < 2) {
      throw("Cannot constrain fit to dimension: No or too few finite estimates available: ", length(x));
    }

    # To find f^{-1}() s.t. x1 = f(x0), we fit x0 = f^{-1}(x1) on (x1,x0).
    fitT <- smooth.spline(lambda, x, ...);
    rm(x,lambda);
    verbose && exit(verbose);

    verbose && enter(verbose, "Adjusting principal-curve fit accordingly");
    for (cc in seq(length=p)) {
      x <- fit$s[,cc, drop=TRUE];
      ok <- whichVector(is.finite(x));
      x[ok] <- predict(fitT, x=x[ok])$y;
      fit$s[,cc] <- values;
    }
    attr(fit, "fitT") <- fitT;
    rm(fitT);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  if (anyMissing) {
    values <- matrix(as.double(NA), nrow=n, ncol=p);
    values[keep,] <- fit$s;
    fit$s <- values;
    values <- rep(as.double(NA), times=n);
    for (ff in c("tag", "lambda")) {
      values[keep] <- fit[[ff]];
      fit[[ff]] <- values;
    }
  }

  verbose && exit(verbose);

  class(fit) <- c("PrincipalCurve", class(fit));

  fit;
}) # fitPrincipalCurve()



###########################################################################
# HISTORY:
# 2008-10-08
# o Now the fitted object is of class PrincipalCurve that extends the
#   princurve::principal.curve class.
# 2008-10-07
# o Added Rdoc comments and an example.
# o Removed implementation for data.frame:s.
# 2008-10-03
# o Added argument 'fixDimension'.
# 2008-05-27
# o Added fitPrincipalCurve().
# o Created.
###########################################################################
