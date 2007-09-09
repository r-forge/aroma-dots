library(aroma.affymetrix)

# source("../aroma.affymetrix/R/CopyNumberSegmentationModel.R");
# source("../aroma.affymetrix/R/GladModel.R");
# source("../aroma.affymetrix/R/CbsModel.R");
# source("../aroma.affymetrix/R/profileCGH.drawCnRegions.R");
# source("../aroma.affymetrix/R/DNAcopy.drawCnRegions.R");
# source("../aroma.affymetrix/R/CopyNumberRegions.R");
# source("../aroma.affymetrix/R/RawCopyNumbers.R");

log <- Verbose(threshold=-4, timestamp=TRUE);

dataSetName <- "Affymetrix_2004-100k_trios,testSet";
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
  plm <- RmaCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, 
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
# Set up list of CN segmenation models
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnsList <- list(
  "glad" = GladModel(cesFlnList),
   "cbs" = CbsModel(cesFlnList)
);
print(cnsList);

# Fit a subset
lapply(cnsList, FUN=fit, arrays=1, chromosomes=19, verbose=log);

# There is a CN deletion on chr 2 @ 83.0Mb in NA06985.
lapply(cnsList, FUN=function(cns) {
  ce <- ChromosomeExplorer(cns);
  process(ce, arrays=1, chromosomes=c(2,19:23), verbose=log);
})

# There is a CN deletion on chr 22 @ 21Mb in NA06994.
lapply(cnsList, FUN=function(cns) {
  ce <- ChromosomeExplorer(cns);
  process(ce, arrays=4, chromosomes=c(2,19:23), verbose=log);
})

