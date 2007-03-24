library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-5);
timestampOn(verbose);

name <- "Affymetrix_2006-HapMap270";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

for (chipType in chipTypes) {
  cs0 <- AffymetrixCelSet$fromName(name, chipType=chipType);

  # Get all parents in the CEU population
  population <- sapply(cs0, FUN=getAttribute, "population");
  fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherid"));
  idxs <- which((population == "CEU") & is.na(fatherID));

  cs <- extract(cs0, idxs);
  print(cs);

  qn <- QuantileNormalization(cs);
  csN <- process(qn, verbose=verbose);
  print(csN);

  plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
  fit(plm, verbose=verbose);
  ces <- getChipEffectSet(plm);
  print(ces);

  fln <- FragmentLengthNormalization(ces);
  cesN <- process(fln, verbose=verbose);
  print(cesN);
}
