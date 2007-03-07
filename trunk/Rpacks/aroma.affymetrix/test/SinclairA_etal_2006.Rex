#savehistory();
library(aroma.affymetrix);
source("init.R");
#source("ui.R");

# Main data set
dataSetName <- "AGRF_2007a";

# Use a subset of the samples
patients <- c(sprintf("%02d", 1:10), sprintf("A%d", 1:5));
controls <- c("C1", "C3", "C7", "S1", "C358,control(|,.*)");
samples <- c(patients, controls);

# Change name of analysis data set
dataSetAlias <- "SinclairA_etal_2006";


chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  # Get main data set
  path <- file.path("rawData", dataSetName, chipType);
  ds <- AffymetrixCelSet$fromFiles(path);

  # Remove duplicates
  dups <- isDuplicated(ds);
  if (any(dups))
    ds <- extract(ds, !dups);

  # Get data set to analyze
  idxs <- indexOf(ds, samples);
  ds <- extract(ds, na.omit(idxs));
  setAlias(ds, dataSetAlias);

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

#  rs <- calculateResiduals(plm, verbose=verbose);
#  ae <- ArrayExplorer(rs);
#  setColorMaps(ae, "log2,rainbow");
#  process(ae, verbose=verbose);

  fln <- FragmentLengthNormalization(ces);
  print(fln);
  print(log, fln);
  cesN <- process(fln, verbose=log);
  print(cesN);
  print(log, cesN);

  cesList[[chipType]] <- cesN;
} # for (chipType ...)

glad <- GladModel(cesList);

ce <- ChromosomeExplorer(glad);
setArrays(ce, patients);
print(ce);

logFile <- sprintf("%s.log", getName(ce));
log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
timestampOn(log);

print(log, ce);
process(ce, verbose=log);
