library(aroma.affymetrix)
source("init.R")

chipTypes <- c("BI_SNP");
dataPath <- "normQuantile";

if (!exists("models", mode="list"))
  models <- list();

for (chipType in chipTypes) {
  model <- models[[chipType]];
  if (!is.null(model))
    next;

  verbose && enter(verbose, "Define dataset");

  # Test set
  dataSetName <- "Affymetrix_2006b-500k_v2,a";
  path <- filePath(dataPath, dataSetName, chipType);
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  verbose && exit(verbose);

  # Define model
  verbose && enter(verbose, "Define model");
  model <- RmaCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);

  verbose && enter(verbose, "Pre-creates probe affinities");
  paf <- getProbeAffinities(model);
  verbose && exit(verbose);

  verbose && enter(verbose, "Pre-creates chip effects");
  ces <- getChipEffects(model);
  verbose && exit(verbose);

  verbose && exit(verbose);

  models[[chipType]] <- model;
 
  # Not needed anymore
  rm(path, dataSetName, ds, dsR, model, ces, paf);
} # for (chipType ...)
gc()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit copy number model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=models)) {
  model <- models[[kk]];
  ds <- getDataSet(model);
  verbose && enter(verbose, "Data set: ", getName(ds));
  
  uDone <- fit(model, moreUnits=30, verbose=verbose);
  nDone <- length(uDone);

  # Update average chip effects
  if (nDone > 0) {
    ces <- getChipEffects(model, verbose=verbose);
    ceAvg <- getAverageFile(ces, units=uDone, force=TRUE, verbose=less(verbose));
    ceAvg$mergeStrands <- model$mergeStrands;
    ceAvg$combineAlleles <- model$combineAlleles;
  }

  verbose && printf(verbose, "Fitted %d SNPs\n", nDone);

  gc();
  verbose && exit(verbose);
} # for (kk ...)


