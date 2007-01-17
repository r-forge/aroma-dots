savehistory();
library(aroma.affymetrix);

dataSetName <- "SinclairA_etal_2006";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  path <- file.path("rawData", dataSetName, chipType);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  # Create a log file
  logFile <- sprintf("%s,%s.log", getName(ds), chipType);
  log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
  timestampOn(log);

  cdf <- getCdf(ds);
  print(cdf);

  gi <- getGenomeInformation(cdf);
  print(gi);

  qn <- QuantileNormalization(ds);
  print(qn);
  dsN <- process(qn, verbose=log);
  print(dsN);

  plm <- RmaCnPlm(dsN, mergeStrands=TRUE, combineAlleles=TRUE);
  print(plm);
  fit(plm, verbose=log);
  ces <- getChipEffects(plm);
  print(ces);

  fln <- FragmentLengthNormalization(ces);
  print(fln);
  cesN <- process(fln, verbose=log);
  print(cesN);

  cesList[[chipType]] <- cesN;
} # for (chipType ...)


# Create a log file
logFile <- sprintf("%s.log", dataSetName);
log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
timestampOn(log);

glad <- GladModel(cesList);
print(glad);
ce <- ChromosomeExplorer(glad);
print(ce);

# Analyze a subset of the data
process(ce, arrays=8, chromosomes=c(7,20:22), verbose=log);

# Whole genome analysis - takes many many hours!
process(ce, arrays=8, chromosomes=c(7,20:22), verbose=log);
