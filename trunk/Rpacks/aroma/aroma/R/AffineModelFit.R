#########################################################################/**
# @RdocClass AffineModelFit
#
# @title "Model fit of affine microarray models"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{a,b,adiag}{Parameters.}
#  \item{residuals}{.}
#  \item{eigen}{.}
#  \item{converged}{If @TRUE, it indicates that the estimate converged. Otherwise, it did not.}
#  \item{nbrOfIterations}{Number of iterations until convergence or stop.}
#  \item{t,t0}{.}
#  \item{fitted}{a @vector of fitted values.}
#  \item{weights}{(only for weighted fits) a @vector the weights used for 
#    the fit. If a iterative re-weighted method is used, these are the 
#    weights used in the last iteration.}
#  \item{y}{if requested, the response @matrix used.}
#  \item{...}{Other named arguments stored.}
# }
#
# \value{
#   Returns a @list containing the fitted model.
# }
#
# @author
#
# \references{
# }
#
# \seealso{
#   @see "RGData.calibrateMultiscan"
#   @see "RGData.normalizeAffine"
#   @seeclass
# }
#*/######################################################################### 
setConstructorS3("AffineModelFit", function(a, b, adiag, fitted=NULL, weights=NULL, residuals=NULL, eigen=NULL, y=NULL, converged=NULL, nbrOfIterations=NULL, t0=NULL, t=NULL, ...) {
  if (missing(a)) {
    fit <- list()
    class(fit) <- "AffineModelFit";
    return(fit);
  }

  if (length(a) != length(b)) {
     throw(sprintf("The length of parameter vector 'a' (%d) and 'b' (%d) are not equal.", length(a), length(b)));
  }
  if (length(a) != length(adiag)) {
     throw(sprintf("The length of parameter vector 'a' (%d) and 'adiag' (%d) are not equal.", length(a), length(adiag)));
  }

  if (any(!isZero(diff(adiag)))) {
     throw("The elements of parameter vector 'adiag' are not all equal: ", paste(adiag, collapse=", "));
  }

  # The number of scans
  nbrOfScans <- length(a);

  N <- NULL;
  if (!is.null(fitted)) {
    if (!is.vector(fitted))
      throw("Argument 'fitted' is not a vector.");
    N <- length(fitted);
  }

  if (!is.null(weights)) {
    if (!is.vector(weights))
      throw("Argument 'weights' is not a vector.");
    if (!is.null(N) && N != length(weights))
      throw("Expected arguments 'weights' to be a vector of length ", N, ": ", length(weights));
    N <- length(weights);
  }

  if (!is.null(residuals)) {
    if (!is.vector(residuals))
      throw("Argument 'residuals' is not a vector.");
    if (!is.null(N) && N != length(residuals))
      throw("Expected arguments 'residuals' to be a vector of length ", N, ": ", length(residuals));
    N <- length(residuals);
  }

  if (!is.null(eigen)) {
    if (!is.matrix(eigen) || diff(dim(eigen)) != 0)
      throw("Argument 'eigen' is not a square matrix.");
    if (nrow(eigen) != nbrOfScans)
      throw(sprintf("Argument 'eigen' is not a %dx%d matrix: %d", nbrOfScans, nbrOfScans, nrow(eigen)));
  }

  if (!is.null(y)) {
    if (!is.matrix(y))
      throw("The response variable 'y' is not a matrix.");
    if (ncol(y) != nbrOfScans) {
      throw(sprintf("The number of scans (%d) does not match the number of columns (%d) in the response matrix.", nbrOfScans, ncol(y)));
    }
  }
 
  # The shortest distance between the fitted line and the diagonal line.
  distance <- sqrt(sum((a - adiag)^2));

  names(a) <- NULL;
  names(b) <- NULL;
  names(adiag) <- NULL;
  coefficients <- c(a=a, b=b, adiag=adiag);

  fit <- list(
    coefficients    = coefficients,
    fitted          = fitted,
    residuals       = residuals,
    weights         = weights,
    eigen           = eigen,
    converged       = converged, 
    nbrOfScans      = nbrOfScans, 
    nbrOfIterations = nbrOfIterations, 
    distance        = distance,
    y               = y,

    # Bootstrap resamples
    t0              = t0,
    t               = t,

    # Elements for backward compatibility
    aPCA            = a,
    k               = b,
    a               = adiag,
    adistance       = distance, 
    x0              = fitted,
    U               = eigen,
    ...
  );
  class(fit) <- "AffineModelFit";
  fit;
})

setMethodS3("as.character", "AffineModelFit", function(this) {
  s <- paste(data.class(this), ":", sep="");
  s <- paste(s, " a=(", paste(format(getA(this)), collapse=","), ")", sep="");
  s <- paste(s, ", b=(", paste(format(getB(this)), collapse=","), ")", sep="");
  s <- paste(s, ", aDiag=(", paste(format(getADiagonal(this)), collapse=","), ")", sep="");
  s <- paste(s, ", aDelta=(", paste(format(getDistance(this)), collapse=","), ")", sep="");
  s;
})

setMethodS3("print", "AffineModelFit", function(x, ...) {
  # To please R CMD check...
  this <- x;

  print(as.character(this), ...);
})

setMethodS3("summary", "AffineModelFit", function(object, ...) {
  # To please R CMD check...
  this <- object;

  summary(this$t, ...);
})

setMethodS3("hasBootstrap", "AffineModelFit", function(this) {
  !is.null(this$t);
})

setMethodS3("fitted", "AffineModelFit", function(object) {
  # To please R CMD check...
  this <- object;

  this$fitted;
})

setMethodS3("residuals", "AffineModelFit", function(object) {
  # To please R CMD check...
  this <- object;

  this$residuals;
})

setMethodS3("getCoefficients", "AffineModelFit", function(this, names=NULL, bootstrap=hasBootstrap(this)) {
  if (is.null(names)) {
    names <- seq(along=this$t0);
  } else {
    K <- this$nbrOfScans;
    names <- paste(rep(names, each=K), 1:K, sep="");
  }

  if (bootstrap) {
    t <- this$t;
    if (is.null(t))
      throw("Can not get bootstrapped coefficients. No bootstrap was performed.");
    t <- t[,names];
    coef <- apply(t, MARGIN=2, FUN=mean);
    sd   <- apply(t, MARGIN=2, FUN=sd);
    attr(coef, "sd*") <- sd;
  } else {
    coef <- this$coefficients;
    coef <- coef[names];
  }
  coef;
})

setMethodS3("coef", "AffineModelFit", function(object, ...) {
  # To please R CMD check...
  this <- object;

  getCoefficients(this, ...);
})


setMethodS3("getA", "AffineModelFit", function(this, ...) {
  getCoefficients(this, names="a", ...);
})

setMethodS3("getADiagonal", "AffineModelFit", function(this, ...) {
  getCoefficients(this, names="adiag", ...);
})

setMethodS3("getB", "AffineModelFit", function(this, ...) {
  getCoefficients(this, names="b", ...);
})


setMethodS3("getDistance", "AffineModelFit", function(this, bootstrap=hasBootstrap(this)) {
  if (bootstrap) {
    t <- this$t;
    if (is.null(t))
      throw("Can not get bootstrapped coefficients. No bootstrap was performed.");
  } else {
    t <- getCoefficients(this);
    t <- t(t);
  }

  K <- this$nbrOfScans;
  a     <- paste("a", 1:K, sep="");
  adiag <- paste("adiag", 1:K, sep="");
  d <- t[,a, drop=FALSE] - t[,adiag, drop=FALSE];
  colnames(d) <- paste("dist", 1:K, sep="");

  d0  <- apply(d, MARGIN=2, FUN=mean);

  if (bootstrap) {
    dsd <- apply(d, MARGIN=2, FUN=sd);
    attr(d0, "sd*") <- dsd;
  }

  d0;
})



############################################################################
# HISTORY:
# 2003-12-29
# o Renamed to AffineModelFit.
# 2003-12-28
# o Created.
############################################################################
