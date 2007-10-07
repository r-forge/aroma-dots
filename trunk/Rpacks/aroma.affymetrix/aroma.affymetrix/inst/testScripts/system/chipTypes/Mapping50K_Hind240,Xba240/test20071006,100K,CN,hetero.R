library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
#chipTypes <- chipTypes[2];

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRawList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
  stopifnot(identical(getNames(cs), sampleNames));
  csRawList[[chipType]] <- cs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csRawList;
csAccList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];

  # Identify cells *not* for chromosome X
  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);
  cells <- getCellIndices(cdf, units=units, useNames=FALSE, unlist=TRUE);
  cells <- setdiff(1:nbrOfCells(cdf), cells);

  acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=cells, tags=c("*", "-X"));
  print(acc);
  csAcc <- process(acc, verbose=log);
  print(csAcc);
  stopifnot(identical(getNames(csAcc), getNames(cs)));
  csAccList[[chipType]] <- csAcc;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csAccList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- HetLogAddCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, 
                                              tags=c("+300", "*", "w"));
  plm$shift <- +300;
  plm$treatNAsAs <- "NA";
  plm$treatNAsAs <- "weighted";
  print(plm);
  fit(plm, ram=1/2, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));
  cesCnList[[chipType]] <- ces;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaList <- list();
cesList <- cesCnList;
for (chipType in names(cesList)) {
  ces <- cesList[[chipType]];
  theta <- extractMatrix(ces, verbose=log);
  print(summary(theta));
  thetaList[[chipType]] <- theta;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesCnList <- cesCnList;
cesFlnList <- list();
for (chipType in names(csList)) {
  ces <- cesCnList[[chipType]];
  fln <- FragmentLengthNormalization(ces);
#  excludeChrXFromFit(fln);  # TO DO
  print(fln);
  cesFln <- process(fln, verbose=verbose);
  print(cesFln);
  stopifnot(identical(getNames(cesFln), getNames(ces)));
  cesFlnList[[chipType]] <- cesFln;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Glad model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
glad <- GladModel(cesFlnList);
print(glad);

fit(glad, arrays=1, chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(glad);
print(ce);
process(ce, chromosomes=c(19,22:23), verbose=log);
## process(ce, verbose=log);
