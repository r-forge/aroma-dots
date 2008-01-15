setMethodS3("calibrateMultiscan", "RawData", function(this, ...) {
  rg <- as.RGData(this);
  fit <- calibrateMultiscan(rg, ...);
  for (ch in c("R", "G"))
    this[[ch]] <- rg[[ch]];
  invisible(fit);
}) # calibrateMultiscan()



setMethodS3("normalizeAffine", "RawData", function(this, ...) {
  rg <- as.RGData(this);
  fit <- normalizeAffine(rg, ...);
  for (ch in c("R", "G"))
    this[[ch]] <- rg[[ch]];
  invisible(fit);
}) # normalizeAffine()



setMethodS3("normalizeWithinSlide", "RawData", function(this, method, fields=c("R", "G", "Rb", "Gb"), lowess=NULL, ...) {
  # Argument 'method':
  if (missing(method))
    throw("Argument 'method' is missing.");

  methods <- rep(method, length.out=nbrOfSlides(this));

  if (any(!is.element(method, c("printorder"))))
    throw("normalizeWithinSlide.RawData: Unknown normalization method: \"", method, "\".");

  # Verify that lowess is a line with fields $x and $y
  if (length(lowess) > 0) {
    if ( any(!is.element(c("x","y"), names(lowess))) ) {
      throw("normalizeWithinSlide.RawData: Can not perform global lowess normaliziation: argument 'lowess' does not contain the field 'x' or 'y'.");
    } else if (any(!is.element(methods, c("printorder")))) {
      throw("normalizeWithinSlide.RawData: When argument 'lowess' is specified, the method must be 'printorder'.");
    }
  }
  
  # Asserts that all needed parameters are specified...
  layout <- getLayout(this);
  if (any(is.element(methods, c("printorder"))) && is.null(layout))
    throw("normalizationWithinSlide.RawData: Can not normalize. Layout object is missing.");

  K <- length(lowess);
  if (K == 1 && all(is.element(c("x","y"), names(lowess))))
    lowess <- list(lowess);

  # Normalize slide by slide.
  for (kk in 1:nbrOfSlides(this)) {
    if (methods[kk] == "printorder") {
      tidx <- toPrintorderMatrix(layout);
      line <- if (length(lowess) > 0) lowess[[(kk-1) %% K + 1]] else NULL;
      description <- "printorder";
      # Order the data in print order
      for (field in fields) {
        X <- this[[field]][tidx];
        idx <- 1:length(tidx);
        # Now, do regular lowess normalization
        ind <- is.na(X) | is.infinite(X);
        if (is.null(line))
          line <- lowess(idx[!ind], X[!ind], ...);
        X[!ind] <- X[!ind] - approx(line, xout=idx[!ind], ties=mean)$y;
        # Update the normalized values
        this[[field]][tidx] <- X;
      } # for (field in ...)
    }
  } # for (kk in 1:nbrOfSlides(this))

  clearCache(this); 

  invisible(this);
}) # normalizeWithinSlide()



setMethodS3("normalizeGenewise", "RawData", function(this, fields=c("R", "G", "Rb", "Gb"), bias=c(10,10), scale=1, ...) {
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})


setMethodS3("normalizeQuantile", "RawData", function(this, 
                                                  fields=c("R", "G"), ...) {
  normalizeQuantile(this, fields=fields, ...);
}) # normalizeQuantile()






############################################################################
# HISTORY:
# 2006-02-08
# o The code for normalizeWithinSlide() assumed 'method' was a single value,
#   but it may be a vector which then generates a warning.
# 2005-03-23
# o Updated normalizeWithinSlide() so that approx() does not give warnings
#   about 'Collapsing to unique x values' when doing lowess normalization.
# 2004-06-27
# o Added calibrateMultiscan() and normalizeAffine().
# o Extracted from RawData.R. 
# 2002-10-24
# o Added normalizeQuantile().
# 2002-04-21
# o Added trial version of normalizeGenewise()
# 2001-08-09
# o Added normalizeWithinSlide.
############################################################################

