setMethodS3("normalizePrincipalCurve", "matrix", function(x, ..., returnFit=FALSE) {
  # Fit principal curve
  fit <- fitPrincipalCurve(x, ...);

  # Flip direction of 'lambda'?
  rho <- cor(fit$lambda, x[,1], use="complete.obs");
  flip <- (rho < 0);

  dx <- (fit$s - x);
  xN <- fit$lambda + dx;
  if (flip) {
    xN <- -xN;
  }

  # Same center for each column
  for (cc in seq(length=ncol(x))) {
    mu <- median(x[,cc], na.rm=TRUE);
    muN <- median(xN[,cc], na.rm=TRUE);
    xN[,cc] <- xN[,cc] - (muN-mu);
  }

  # Return fit?
  if (returnFit)
    attr(xN, "fit") <- fit;

  xN;
}) # normalizePrincipalCurve()


setMethodS3("normalizePrincipalCurve", "data.frame", function(x, ...) {
  xN <- normalizePrincipalCurve(as.matrix(x), ...);
  xN <- as.data.frame(xN);
})



setMethodS3("makeSmoothSplinePredict", "numeric", function(x, y, df=5, ...) {
  # Argument 'y':
  if (length(x) != length(y)) {
    throw("Argument 'y' is of a different length than 'x'");
  }

  # Identify finite (x,y) pairs
  ok <- which(is.finite(x) & is.finite(y));
  x <- x[ok];
  y <- y[ok];
  rm(ok);

  specs <- list(
    xRange = range(x),
    yRange = range(y)
  );

  # Fit smooth function
  fit <- smooth.spline(x,y, df=df, ...);
  rm(x,y);

  # Create predict() function that handles missing values.
  predFcn <- function(x, ...) {
    yPred <- rep(as.double(NA), length(x));
    ok <- which(!is.na(x));  # Allows for -/+Inf:s though
    yPred[ok] <- predict(fit, x=x[ok], ...)$y;
    yPred;
  } # predFcn()

  attr(predFcn, "fit") <- fit;
  attr(predFcn, "specs") <- specs;

  predFcn;
}) # makeSmoothSplinePredict()



###########################################################################
# HISTORY:
# 2008-05-27
# o Added normalizePrincipalCurve().
# o Created.  Will probably end up in aroma.light.
###########################################################################
