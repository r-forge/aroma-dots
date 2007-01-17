savehistory();
library(aroma.affymetrix);
source("init.R");
source("ui.R");

paths <- c("rawData");
pattern <- "_(199|200)[0-9]";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}

cesList <- list();
for (ds in dataSets) {
  # Remove duplicates
  dups <- isDuplicated(ds);
  if (any(dups))
    ds <- extract(ds, !dups);
  print(ds);

  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);

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
