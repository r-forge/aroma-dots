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
    fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherID"));
    idxs <- which((population == "CEU") & is.na(fatherID));
  
    cs <- extract(cs0, idxs);
    setAlias(cs, sprintf("%s.CEU.founders", getName(cs0)));
    rm(cs0);
    n23 <- as.integer(sapply(cs, FUN=getAttribute, "n23"));
  }
  print(cs);

  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);

  subsetToFit <- seq_len(nbrOfUnits(cdf));
  subsetToFit <- setdiff(subsetToFit, units);

#  csC <- process(acc, verbose=verbose);
#  bg <- RmaBackgroundCorrection(cs);

  acc <- AllelicCrosstalkCalibration(cs);
  print(acc);
  csC <- process(acc, verbose=verbose);
  print(csC);

#  qn <- DChipQuantileNormalization(csC, typesToUpdate="pm", 
#                                                 tags=c("*", "pm", "-X"));
#  excludeChrXFromFit(qn);
#  print(qn);
#  csN <- process(qn, verbose=verbose);
#  print(csN);

  gc();

#  plm <- MbeiCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("-100","*"));
  plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("+300","*"));
  plm$treatNAsAs <- "NA";
  plm$shift <- +300;
  print(plm);
  ces <- getChipEffectSet(plm);
  fit(plm, ram=2.5, verbose=verbose);

#  fit(plm, units=units, ram=2.5, verbose=verbose);

  fln <- FragmentLengthNormalization(ces, subsetToFit=subsetToFit, tags=c("*", "-X"));
  print(fln);
  cesN <- process(fln, verbose=verbose);
  print(cesN);

#  gcn <- GcContentNormalization(cesN, subsetToFit=subsetToFit, tags=c("*", "-X"));
#  cesN2 <- process(gcn, verbose=verbose);

#  ces <- cesN;

  qn2 <- DChipQuantileNormalization(cesN, tags=c("*", "-X"));
  excludeChrXFromFit(qn2);
  print(qn2);

  cesQN <- process(qn2, verbose=verbose);
#  print(csN);

stop()

  print(ces);
  if (TRUE) {
    rm(theta); 
    theta <- extractMatrix(ces, units=units, verbose=verbose);
    saveObject(list(theta=theta), sprintf("%s.RData", paste(getTags(ces), collapse=",")), compress=FALSE);
    rm(theta); 
  }

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
