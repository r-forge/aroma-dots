library(aroma.affymetrix)

pngDev <- System$findGraphicsDevice();

chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

dataPath <- "normQuantile";

if (!exists("models", mode="list"))
  models <- list();

for (chipType in chipTypes) {
  model <- models[[chipType]];
  if (!is.null(model))
    next;

  verbose && enter(verbose, "Define dataset");

  # Test set
  dataSetName <- "SlaterH_etal_2004,QN";
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
allChromosomes <- c(1:22, "X");
#allChromosomes <- c(11, 14, 21);
fitModels <- TRUE;
#fitModels <- FALSE;
if (fitModels) {
  for (kk in seq(along=models)) {
    model <- models[[kk]];
    ds <- getDataSet(model);
    verbose && enter(verbose, "Data set: ", getName(ds));
    
    for (chr in allChromosomes) {
      verbose && enter(verbose, "Chromosome ", chr);
  
      verbose && enter(verbose, "Retrieving genome information");
      cdf <- getCdf(model);
      gi <- getGenomeInformation(cdf);
      units <- getUnitIndices(gi, chromosome=chr);
      nunits <- length(units);
      verbose && exit(verbose);
    
      verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
      uDone <- fit(model, units=units, moreUnits=30, verbose=verbose);
      nDone <- length(uDone);
    
      # Update average chip effects
      if (nDone > 0) {
        ces <- getChipEffects(model, verbose=verbose);
        ceAvg <- getAverageFile(ces, units=uDone, force=TRUE, verbose=less(verbose));
        ceAvg$mergeStrands <- model$mergeStrands;
        ceAvg$combineAlleles <- model$combineAlleles;
      }
    
      verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
      verbose && exit(verbose);
      verbose && exit(verbose);
    }

    gc();
    verbose && exit(verbose);
  } # for (chr ...)
} # if (fitModels)


