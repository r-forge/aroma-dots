onParameters <- function(..., pluginVersion="0.3.0") {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Checking plugin version
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (compareVersion(pluginVersion, "0.3.0") < 0)
    throw("Requested plugin version is too old: ", pluginVersion);
  #L-#


  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Pre-processing plugin parameters
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L-#

  # Return updated plugin parameters
  list()
}
