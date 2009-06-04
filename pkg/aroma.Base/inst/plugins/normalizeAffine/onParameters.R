onParameters <- function(method="L1", constraint=0.05, multiassay=c("singleassay", "multiassay", "multiassayPerChannel", "multiassayWithReferenceChannel"), referenceChannel=c("1","2"), ..., pluginVersion="0.5") {
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
  validValues <- eval(formals(normalizeAffine.matrix)[["method"]]);
  method <- match.arg(method, validValues);

  # Argument 'bandwidth':
  constraint <- as.double(constraint);
  if (is.na(constraint))
    throw("Parameter 'constraint' is NA.");
  if (constraint <= 0)
    constraint <- 0.05;

  # Argument 'multiassay':
  multiassay <- match.arg(multiassay);

  # Argument 'referenceChannel':
  referenceChannel <- unlist(strsplit(referenceChannel, split="[\t, ]"));
  referenceChannel <- as.integer(referenceChannel);
  if (any(is.na(referenceChannel))) {
    throw("Parameter 'referenceChannel' is NA: ", 
                                    paste(referenceChannel, collapse=", "));
  }
  if (any(referenceChannel < 1 | referenceChannel > 2)) {
    throw("Parameter 'referenceChannel' is out of range: ", 
                                    paste(referenceChannel, collapse=", "));
  }

  #L-#

  # Return updated plugin parameters
  list(
    method=method, 
    constraint=constraint, 
    multiassay=multiassay, 
    referenceChannel=referenceChannel
  );
}
