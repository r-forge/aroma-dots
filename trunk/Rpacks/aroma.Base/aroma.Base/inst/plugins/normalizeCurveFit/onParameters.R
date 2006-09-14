onParameters <- function(method="lowess", bandwidth=NULL, ..., pluginVersion="1.0") {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Checking plugin version
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (compareVersion(pluginVersion, "0.5") > 0)
    throw("Requested plugin version is not available: ", pluginVersion);

  if (compareVersion(pluginVersion, "0.5") < 0)
    throw("Requested plugin version is too old: ", pluginVersion);
  #L-#

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Pre-processing plugin parameters
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'methods':
  validValues <- eval(formals(normalizeCurveFit.matrix)[["method"]]);
  method <- match.arg(method, validValues);

  # Argument 'bandwidth':
  bandwidth <- as.double(bandwidth);
  if (is.na(bandwidth) || bandwidth <= 0)
    bandwidth <- NULL;
  #L-#

  # Return updated plugin parameters
  list(
    method=method,
    bandwidth=bandwidth
  );
}
