library(aroma.affymetrix);
source("../test/init.R");
verbose <- Arguments$getVerbose(-5);
timestampOn(verbose);


name <- "Affymetrix_2006-HapMap270";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  path <- filePath("rawData", name, chipType, expandLinks="any");
  cs0 <- AffymetrixCelSet$fromFiles(path);

  # Get all parents in the CEU population
  population <- sapply(cs0, FUN=getAttribute, "population");
  fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherid"));
  idxs <- which((population == "CEU") & is.na(fatherID));

  cs <- extract(cs0, idxs);
  setAlias(cs, sprintf("%s.CEU.founders", getName(cs0)));
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

  cesList[[chipType]] <- cesN;
}

# Set sample annotation attributes
saf <- SampleAnnotationSet$fromPath("rawData");
cesList <- lapply(cesList, setAttributesBy, saf);

# Calculate baseline according to Bengtsson H (2007).
refList <- lapply(cesList, calculateBaseline, defaultPloidy=2, verbose=verbose)

glad <- GladModel(cesList, refList);

ce <- ChromosomeExplorer(glad);
process(ce, chromosome=23, verbose=verbose);

