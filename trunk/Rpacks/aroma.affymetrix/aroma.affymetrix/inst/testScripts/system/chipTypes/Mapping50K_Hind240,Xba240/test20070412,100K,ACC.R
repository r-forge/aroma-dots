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
