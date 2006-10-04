pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "normQuantile";

name <- "Affymetrix_2006-CEPH270";
name <- "Affymetrix_2005-500k_data";
name <- "Affymetrix_2004-100k_trios";
#name <- "Affymetrix_2004-100k_trios-calibACT";
#name <- "Affymetrix_2005-500k_data-calibACT";
#name <- "SinclairA_2006-calibACT";
#name <- "SlaterH_etal_2004-calibACT";

if (!exists("dataPath"))
  dataPath <- findCelSet(name=name, paths=paths, minCount=5);
print(dataPath);

if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(dataPath);
  # Get both chiptypes in the chip set
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}
print(raw);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define copy-number PLM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  verbose && enter(verbose, "Defining PLM models");
  models <- vector("list", length(raw));
  names(models) <- names(raw);
  for (kk in seq(along=models)) {
    ds <- raw[[kk]];
    verbose && enter(verbose, "Data set: ", getName(ds));
#    model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
    model <- RmaCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
    models[[kk]] <- model;
  }
  verbose && print(verbose, model);

  verbose && enter(verbose, "Pre-creates probe affinities");
  paf <- getProbeAffinities(model, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Pre-creates chip effects");
  ces <- getChipEffects(model, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && exit(verbose);
}
print(models);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit PLM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

allChromosomes <- c(1:22, "X");
for (chr in allChromosomes) {
  verbose && enter(verbose, "Chromosome ", chr);

  for (kk in seq(along=models)) {
    model <- models[[kk]];
    ds <- getDataSet(model);
    verbose && enter(verbose, "Data set: ", getName(ds));
  
    verbose && enter(verbose, "Retrieving genome information");
    cdf <- getCdf(model);
    gi <- getGenomeInformation(cdf);
    units <- getUnitIndices(gi, chromosome=chr);
    nunits <- length(units);
    verbose && exit(verbose);
  
    verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
    uDone <- fit(model, units=units, moreUnits=1, verbose=verbose);
    nDone <- length(uDone);
  
    verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
    verbose && exit(verbose);
    verbose && exit(verbose);
  }
  
  verbose && exit(verbose);
} # for (chr in c(1:22, "X"))

