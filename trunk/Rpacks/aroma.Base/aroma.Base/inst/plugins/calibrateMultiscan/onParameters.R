onParameters <- function(method="L1", assayNamePattern="^(.+)_(.+)$", groupBy="\\1", reports=c("html", "none"), imageSize=c(420,420), ..., pluginVersion="0.5") {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Checking plugin version
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (compareVersion(pluginVersion, "0.5") < 0)
    throw("Requested plugin version is too old: ", pluginVersion);
  #L-#

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Pre-processing plugin parameters
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'methods':
  validValues <- eval(formals(calibrateMultiscan.matrix)[["method"]]);
  method <- match.arg(method, validValues);

  # Argument 'assayNamePattern':
  assayNamePattern <- Arguments$getCharacter(assayNamePattern);

  # Argument 'groupBy':
  groupBy <- Arguments$getCharacter(groupBy);

  # Argument 'reports':
  reports <- match.arg(reports);

  # Argument 'imageSize':
  imageSize <- Arguments$getIntegers(imageSize, length=2, range=c(16,Inf));
  
  #Lc# Testing pattern 'assayNamePattern': ${assayNamePattern}
  gsub(assayNamePattern, "", "Hello world");

  #Lc# Testing 'groupBy': ${groupBy}
  gsub(assayNamePattern, groupBy, "Hello world");

  #L-#

  # Return updated plugin parameters
  list(
    method=method, 
    assayNamePattern=assayNamePattern,
    groupBy=groupBy,
    reports=reports,
    imageSize=imageSize
  );
}
