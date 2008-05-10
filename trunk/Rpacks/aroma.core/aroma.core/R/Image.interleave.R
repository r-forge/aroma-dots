setMethodS3("interleave", "Image", function(this, what=c("none", "h", "v", "auto"), ...) {
  require("EBImage") || throw("Package not loaded: EBImage.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  safeMeans <- function(x) {
    mean(x[is.finite(x)]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Interleave
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Nothing todo?
  if (what == "none")
    return(this);

  verbose && enter(verbose, "Interleaving image");

  # Get image data
  z <- this@.Data;
  
  # if only PM locations have signal, add a fake row
  if (what == "auto") {
    verbose && enter(verbose, "Infering horizontal, vertical, or no interleaving");
    n <- 2*(nrow(z) %/% 2);
    idxOdd <- seq(from=1, to=n, by=2);
    hOdd <- safeMeans(abs(z[idxOdd,,]));
    hEven <- safeMeans(abs(z[idxOdd+1,,]));
    verbose && printf(verbose, "hOdd=%.2g\n", hOdd);
    verbose && printf(verbose, "hEven=%.2g\n", hEven);
    hRatio <- log(hOdd/hEven);
    verbose && printf(verbose, "hRatio=%.2g\n", hRatio);

    n <- 2*(ncol(z) %/% 2);
#    n <- max(n, 40);  # Infer from the first 40 rows.
    idxOdd <- seq(from=1, to=n, by=2);
    vOdd <- safeMeans(abs(z[,idxOdd,]));
    vEven <- safeMeans(abs(z[,idxOdd+1,]));
    verbose && printf(verbose, "vOdd=%.2g\n", vOdd);
    verbose && printf(verbose, "vEven=%.2g\n", vEven);
    vRatio <- log(vOdd/vEven);
    verbose && printf(verbose, "vRatio=%.2g\n", vRatio);

    what <- "none";
    if (abs(vRatio) > abs(hRatio)) {
      if (abs(vRatio) > 0.25) {
        if (vRatio > 0)
          what <- "v"
        else
          what <- "v";
      }
    } else {
      if (abs(hRatio) > 0.25) {
        if (hRatio > 0)
          what <- "h"
        else
          what <- "h";
      }
    }
    verbose && cat(verbose, "what: ", what);
    verbose && exit(verbose);
  }

  isUpdated <- FALSE;
  if (what == "h") {
    idxOdd <- seq(from=1, to=2*(nrow(z) %/% 2), by=2);
    z[idxOdd,,] <- z[idxOdd+1,];
  } else if (what == "v") {
    idxOdd <- seq(from=1, to=2*(ncol(z) %/% 2), by=2);
    z[,idxOdd,] <- z[,idxOdd+1,];
  } else {
    isUpdated <- FALSE;
  }

  # Update?
  if (isUpdated) {
    this@.Data <- z;
  }

  verbose && exit(verbose);

  this;
})



############################################################################
# HISTORY:
# 2008-05-10
# o BUG FIX: interleave() for Image gave 'Error in z[idxOdd, ] : incorrect 
#   number of dimensions'.  The internal image structure is a 3-dim array.
# 2008-03-14
# o Created from getImage() of AffymetrixCelFile.
############################################################################

