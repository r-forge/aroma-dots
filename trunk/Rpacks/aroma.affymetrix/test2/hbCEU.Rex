library(aroma.affymetrix);
source("../test/init.R");
verbose <- Arguments$getVerbose(-5);
timestampOn(verbose);

name <- "Affymetrix_2006-HapMap270";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

for (chipType in chipTypes[1]) {
  print(chipType);

  if (!exists("cs")) {
    path <- filePath("rawData", name, chipType, expandLinks="any");
    cs0 <- AffymetrixCelSet$fromFiles(path);
  
    # Get all parents in the CEU population
    population <- sapply(cs0, FUN=getAttribute, "population");
    fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherid"));
    idxs <- which((population == "CEU") & is.na(fatherID));
  
    n23 <- as.integer(sapply(cs0, FUN=getAttribute, "n23"));
  
    cs <- extract(cs0, idxs);
    setAlias(cs, sprintf("%s.CEU.founders", getName(cs0)));
    rm(cs0);
  }
  print(cs);

  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);

#  acc <- AllelicCrosstalkCalibration(cs);
#  csC <- process(acc, verbose=verbose);
#  bg <- RmaBackgroundCorrection(cs);
#  print(bg);
#  csC <- process(bg, verbose=verbose);
#  print(csC);

#  qn <- DChipQuantileNormalization(cs);
#  qn <- DChipQuantileNormalization(cs, tags=c("*", "-X"));
#  excludeChrXFromFit(qn);
  qn <- QuantileNormalization(cs);
#  csN2 <- process(qn2, verbose=verbose);
  print(qn);
  csN <- process(qn, verbose=verbose);
  print(csN);

  gc();

  plm <- MbeiCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("*","-100"));
  plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("*","-100"));
  plm <- AffineCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, background=FALSE);
  plm$shift <- -100;
  print(plm);
  ces <- getChipEffectSet(plm);
  fit(plm, units=units, ram=2, verbose=verbose);
  if (FALSE) {
    rm(theta); 
    theta <- extractMatrix(ces, units=units, verbose=verbose);
    saveObject(list(theta=theta), sprintf("%s.RData", getFullName(ces)), compress=FALSE);
  }
  print(ces);

  fln <- FragmentLengthNormalization(ces);
  cesN <- process(fln, verbose=verbose);
  print(cesN);

  if (!plm$mergeStrands) {
    merger <- ChipEffectGroupMerge(cesN);
    setMergeFunction(merger, function(x) { 
      x[1,] <- 2^colMeans(log2(x), na.rm=TRUE); 
      x[-1,] <- 0;
      x;
    })
    cesM <- process(merger, verbose=verbose);
  }

  gc();
}
