setMethodS3("backtransformPrincipalCurve", "matrix", function(X, fit, dimensions=NULL, targetDimension=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'X'
  if (!is.numeric(X)) {
    throw("Argument 'X' is not numeric: ", mode(X));
  }

  dim <- dim(X);
  if (!is.matrix(X)) {
    X <- as.matrix(X);
  }

  # Argument 'fit'
  if (!inherits(fit, "principal.curve")) {
    throw("Argument 'fit' is not a principal.curve object: ", class(fit)[1]);
  }

  # Argument 'dimensions'
  p <- ncol(fit$s);
  if (!is.null(dimensions)) {
    dimensions <- Arguments$getIndices(dimensions, range=c(1, p));
  }

  if (!is.null(targetDimension)) {
    targetDimension <- Arguments$getIndex(targetDimension, range=c(1, p));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Pre-allocate result matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  naValue <- NA;
  mode(naValue) <- mode(X);
  Xhat <- matrix(naValue, nrow=dim[1], ncol=dim[2]);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Transform towards a target dimension?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  hasTargetDimension <- (!is.null(targetDimension));
  if (hasTargetDimension) {
    lambda <- fit$s[,targetDimension];
  } else {
    lambda <- fit$lambda;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Find backtransformations and backtransform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  s <- fit$s;
  if (!is.null(dimensions)) {
    s <- s[,dimensions,drop=FALSE];
  }

    
  for (kk in seq(length=ncol(s))) {
    sKK <- s[,kk];
    fitKK <- smooth.spline(sKK, lambda, ...);

    Xkk <- X[,kk];
    keep <- whichVector(is.finite(Xkk));
    Xkk <- Xkk[keep];
    XhatKK <- predict(fitKK, x=Xkk)$y;
    stopifnot(length(XhatKK) == length(keep));
    Xhat[keep,kk] <- XhatKK;
  }

  rm(sKK, lambda, fitKK, XhatKK, keep, s);

  dim(Xhat) <- dim;
  Xhat;
}) # backtransformPrincipalCurve()


setMethodS3("backtransformPrincipalCurve", "numeric", function(X, ...) {
  X <- as.matrix(X);
  backtransformPrincipalCurve(X, ...);
})

###########################################################################
# HISTORY:
# 2008-10-08
# o Added argument 'targetDimension' to backtransformPrincipalCurve().
# 2008-10-07
# o Created.
###########################################################################
