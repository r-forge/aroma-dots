library(aroma.affymetrix);
#source("../test/init.R");
verbose <- Arguments$getVerbose(-5);
timestampOn(verbose);

name <- "Affymetrix_2006-HapMap270";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

for (chipType in chipTypes[1]) {
  path <- filePath("rawData", name, chipType, expandLinks="any");
  cs0 <- AffymetrixCelSet$fromFiles(path);

  # Get all parents in the CEU population
  population <- sapply(cs0, FUN=getAttribute, "population");
  fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherid"));
  idxs <- which((population == "CEU") & is.na(fatherID));

  n23 <- as.integer(sapply(cs0, FUN=getAttribute, "n23"));

  cs <- extract(cs0, idxs);
  setAlias(cs, sprintf("%s.CEU.founders", getName(cs0)));
  print(cs);
  rm(cs0);

  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);

#  acc <- AllelicCrosstalkCalibration(cs);
#  csC <- process(acc, verbose=verbose);
#  bg <- RmaBackgroundCorrection(cs);
#  csC <- process(bg, verbose=verbose);
#  print(csC);

  qn <- QuantileNormalization(cs);
  csN <- process(qn, verbose=verbose);
  print(csN);

#  plm <- MbeiCnPlm(csN, mergeStrands=FALSE, combineAlleles=TRUE);
  plm <- RmaCnPlm(csN, mergeStrands=FALSE, combineAlleles=TRUE);
  fit(plm, ram=4, verbose=verbose);
  ces <- getChipEffectSet(plm);
  print(ces);

  fln <- FragmentLengthNormalization(ces);
  cesN <- process(fln, verbose=verbose);
  print(cesN);

  merger <- ChipEffectGroupMerge(cesN);
  setMergeFunction(merger, function(x) { 
    x[1,] <- 2^colMeans(log2(x), na.rm=TRUE); 
    x[-1,] <- 0;
    x;
  })
  cesM <- process(merger, verbose=verbose);

  gc();
}
