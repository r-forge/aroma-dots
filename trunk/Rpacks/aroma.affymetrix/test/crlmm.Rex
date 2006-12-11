library(oligo);
library(aroma.affymetrix);
source("init.R");
source("ui.R");

moreUnits <- 30;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "probeData";
pattern <- "^[a-zA-Z0-9_-]*$";
pattern <- NULL;
paths <- "raw";
pattern <- "^[^,]*$";
pattern <- "^[a-zA-Z0-9_-]*$";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  rm(ds);
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");

res <- list()
for (kk in seq(along=dataSets)) {
  ds <- dataSets[[kk]];
  verbose && enter(verbose, sprintf("Data set %d", kk));
  verbose && print(verbose, ds);

  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibration of probe signals for PCR effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sc <- PcrProbeCalibration(ds);
  verbose && print(verbose, sc);
  dsC <- process(sc, verbose=verbose);
  verbose && print(verbose, dsC);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Quantile normalization of PMs toward HapMap references
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hapmapQs <- getReferenceQuantiles(pd);
  norm <- QuantileNormalization(dsC, typesToUpdate="pm", 
         targetDistribution=hapmapQs, tags=c("*", "HapMap"));
  verbose && print(verbose, norm);
  dsN <- process(norm, verbose=verbose);
  verbose && print(verbose, dsN);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit RMA
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  plm <- RmaSnpPlm(dsN, mergeStrands=FALSE);
  verbose && print(verbose, plm);
  units <- fit(plm, moreUnits=moreUnits, verbose=verbose);
  verbose && exit(verbose);

}

verbose && exit(verbose);

