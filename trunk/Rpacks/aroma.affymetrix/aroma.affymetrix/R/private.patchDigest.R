.patchDigest <- function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require("affxparser") || throw("Package not loaded: affxparser");
  env <- as.environment("package:affxparser");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch digest()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If using digest v0.2.3 or before, we need to use the above fix.
  # In digest v0.2.4 this is the default behavior.
  if (compareVersion(packageDescription("digest")$Version, "0.3.0") >= 0) {
    digest.patch <- function(..., skip="auto", ascii=FALSE) {
      digest::digest(..., skip=skip, ascii=ascii);
    }
  } else {
    digest.patch <- .digest.fix023;
  }
  env <- as.environment("package:aroma.affymetrix");
  assign("digest", digest.patch, envir=env);

  invisible();
} # .patchDigest()

############################################################################
# HISTORY:
# 2007-04-04
# o Created.
############################################################################  
