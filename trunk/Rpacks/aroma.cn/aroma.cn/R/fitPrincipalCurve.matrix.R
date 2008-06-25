setMethodS3("fitPrincipalCurve", "matrix", function(x, ..., .useNew=TRUE, verbose=FALSE) {
  require("princurve") || throw("Package not loaded: princurve");

  # The current implementation contains bugs. /HB 2008-05-26
  principal.curve <- principal.curve.hb;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting principal curve");
  n <- nrow(x);
  p <- ncol(x);

  verbose && cat(verbose, "Data size: ", n, "x", p);

  # princurve::principal.curve() does not handle missing values.
  keep <- rep(TRUE, n);
  for (cc in seq(length=p))
    keep <- keep & is.finite(x[,cc]);
  anyMissing <- (!all(keep));
  if (anyMissing)
    x <- x[keep,];

  verbose && cat(verbose, "Data size after removing non-finite data points: ", n-sum(keep), "x", p);


  verbose && enter(verbose, "Calling principal.curve()");
  fit <- principal.curve(x, ...);
  verbose && exit(verbose);

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

  fit;
}) # fitPrincipalCurve()


setMethodS3("fitPrincipalCurve", "data.frame", function(x, ...) {
  fitPrincipalCurve(as.matrix(x), ...);
})


###########################################################################
# HISTORY:
# 2008-05-27
# o Added fitPrincipalCurve().
# o Created.
###########################################################################
