setMethodS3("normalizeFragmentLength", "default", function(y, fragmentLengths, targetFcn=NULL, subsetToFit=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  y <- Arguments$getDoubles(y, disallow=NULL);
  
  # Argument 'fragmentLengths':
  fragmentLengths <- Arguments$getDoubles(fragmentLengths, length=length(y), disallow=NULL);

  # Argument 'targetFcn':
  if (!is.null(targetFcn)) {
    if (!is.function(targetFcn)) {
      throw("Argument 'targetFcn' must be a function or NULL: ", class(targetFcn)[1]);
    }
  }
  
  # Argument 'subsetToFit':
  if (!is.null(subsetToFit)) {
    subsetToFit <- Arguments$getIndices(subsetToFit, length=length(y));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate normalization function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit smooth curve
  ok <- (is.finite(fragmentLengths) & is.finite(y));
  if (!is.null(subsetToFit)) {
    ok[-subsetToFit] <- FALSE;
  }
  suppressWarnings({
    fit <- lowess(fragmentLengths[ok], y[ok], ...);
    class(fit) <- "lowess";
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dy <- predict(fit, newdata=fragmentLengths);

  # Normalize toward a target function?
  if (!is.null(targetFcn)) {
    yTargetPred <- targetFcn(fragmentLengths);
    dy <- dy - yTargetPred;
  }

  # Transform signals
  ok <- is.finite(dy);
  y[ok] <- y[ok] - dy[ok];
  
  y;
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-11-28
# o Created.
############################################################################
