setMethodS3("backtransformPrincipalCurve", "matrix", function(X, fit, dimensions=NULL, ...) {
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
  if (!is.null(dimensions)) {
    dimensions <- Arguments$getIndices(dimensions, range=c(1, dim[2]));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Pre-allocate result matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  naValue <- NA;
  mode(naValue) <- mode(X);
  Xhat <- matrix(naValue, nrow=dim[1], ncol=dim[2]);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Find backtransformations and backtransform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  s <- fit$s;
  if (!is.null(dimensions)) {
    s <- s[,dimensions,drop=FALSE];
  }


  for (kk in seq(length=ncol(s))) {
    sKK <- s[,kk];
    lambdaKK <- fit$lambda;
    fitKK <- smooth.spline(sKK, lambdaKK, ...);
    XhatKK <- predict(fitKK, x=X[,kk])$y;
    stopifnot(length(XhatKK) == dim[1]);
    Xhat[,kk] <- XhatKK;
  }
  rm(sKK, lambdaKK, fitKK, XhatKK, s);

  dim(Xhat) <- dim;
  Xhat;
}) # backtransformPrincipalCurve()



###########################################################################
# HISTORY:
# 2008-10-07
# o Created.
###########################################################################
