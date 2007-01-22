savehistory();
library(aroma.affymetrix);
source("init.R");
source("ui.R");

dataSetName <- "AGRF_2007a-Anonymous";

chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  path <- file.path("rawData", dataSetName, chipType);
  ds <- AffymetrixCelSet$fromFiles(path);

  # Remove duplicates
  dups <- isDuplicated(ds);
  if (any(dups))
    ds <- extract(ds, !dups);
  print(ds);

  # Create a log file
  logFile <- sprintf("%s,%s.log", getName(ds), chipType);
  log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
  timestampOn(log);
  print(log, ds);

  qn <- QuantileNormalization(ds);
  print(qn);
  print(log, qn);

  dsN <- process(qn, verbose=log);
  print(dsN);
  print(log, dsN);

  plm <- RmaCnPlm(dsN, mergeStrands=TRUE, combineAlleles=TRUE);
  print(plm);
  print(log, plm);
  fit(plm, verbose=log);
  ces <- getChipEffects(plm);
  print(ces);
  print(log, ces);

  fln <- FragmentLengthNormalization(ces);
  print(fln);
  print(log, fln);
  cesN <- process(fln, verbose=log);
  print(cesN);
  print(log, cesN);

  cesList[[chipType]] <- cesN;
} # for (chipType ...)
