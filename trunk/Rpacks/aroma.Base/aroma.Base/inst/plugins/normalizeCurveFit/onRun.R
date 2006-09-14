onRun <- function(data, method="lowess", bandwidth=NULL, ..., pluginVersion="1.0") {
  spots <- getSections(data, "spots");

  dprogress <- 95/length(spots);

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #L+# Normalizing ${length(spots)} assays
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (spot in spots) {
    #L+# Normalizing assay '${getAssays(spot)[1]}'

    #Lc# Extracting intensities
    X <- getData(spot, fields=c("intensity1", "intensity2"));
    X <- as.matrix(X);

    #Lc# Calling normalization function
    X <- normalizeCurveFit(X, method=method, bandwidth=bandwidth);
    colnames(X) <- c("intensity1", "intensity2");

    #Lc# Updating intensities
    setDataFields(spot, values=X);

    increase(progress, dprogress);
    #L-#
  }
  #L-#

  # Returning modified data
  spots;
}
