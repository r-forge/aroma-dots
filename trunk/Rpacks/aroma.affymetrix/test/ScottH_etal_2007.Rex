savehistory();
library(aroma.affymetrix);
source("init.R");

dataSetName <- "AGRF_2007a-Anonymous";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  path <- file.path("plmData", dataSetName, chipType);
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


# Create a log file
glad <- GladModel(cesList);
print(glad);

ce <- ChromosomeExplorer(glad);
setAlias(ce, "ScottH_etal_2007");
print(ce);

logFile <- sprintf("%s.log", getName(ce));
log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
timestampOn(log);
print(log, glad);
print(log, ce);

arrays <- c("1-288", "2-437", "3-575");
arrays <- which(getArrays(glad) %in% arrays);
process(ce, arrays=arrays, verbose=log);
