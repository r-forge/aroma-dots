savehistory();
library(aroma.affymetrix);
source("init.R");

dataSetName <- "SinclairA_etal_2006";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

#dataSetName <- "SlaterH_etal_2004";
#chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

cesList <- list();
for (chipType in chipTypes) {
  # Create a log file
  logFile <- sprintf("%s,%s.log", dataSetName, chipType);
  log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
  timestampOn(log);

  path <- file.path("rawData", dataSetName, chipType);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

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

process(ce, arrays=8, chromosomes=c(7,20:22), verbose=log);
