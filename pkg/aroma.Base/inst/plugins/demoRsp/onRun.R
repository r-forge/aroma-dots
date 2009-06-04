onRun <- function(data, ..., pluginVersion="0.5") {
  assays <- getSections(data, "spots");
  nbrOfAssays <- length(assays);

  #Lc# pluginPath=${pluginPath}

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Listing assays
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (assay in assays) {
  }
  gc();

  #L-#

  # Returning nothing.
  NULL;
}
