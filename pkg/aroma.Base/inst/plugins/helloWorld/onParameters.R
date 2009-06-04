onParameters <- function(welcomeMessage="Hello world!", ..., pluginVersion="0.3.0") {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Checking plugin version
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (compareVersion(pluginVersion, "0.3.0") < 0)
    throw("Requested plugin version is too old: ", pluginVersion);
  #L-#


  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Pre-processing plugin parameters
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nchar(welcomeMessage) == 0) {
    #Lc# Empty welcome message detected. Using default message.
    welcomeMessage <- formals(sys.function())[["welcomeMessage"]];
  }

  if (nchar(welcomeMessage) > 256) {
    throw("Ridiculous long welcome message: ", welcomeMessage);
  }
  #L-#

  # Return updated plugin parameters
  list(welcomeMessage=welcomeMessage)
}
