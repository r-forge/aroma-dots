.patchAffxparser <- function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require("affxparser") || throw("Package not loaded: affxparser");
  env <- as.environment("package:affxparser");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch findCdf()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  environment(findCdf.patch) <- env;
  unlockBinding("findCdf", env);
  assignInNamespace("findCdf", findCdf.patch, ns="affxparser", envir=env);
  assign("findCdf", findCdf.patch, envir=env);
  lockBinding("findCdf", env);

  invisible();
} # .patchAffxparser()

############################################################################
# HISTORY:
# 2007-02-12
# o Created.
############################################################################  
