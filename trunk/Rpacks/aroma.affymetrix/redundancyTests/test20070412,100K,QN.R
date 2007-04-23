library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

data - 57448.854736
dataSetName <- "Affymetrix_2004-100k_trios,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
chipTypes <- chipTypes[1];

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
# Quantile normalization tests #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csQnList <- list();
csList <- csRawList;
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  qn <- QuantileNormalization(cs);
  print(qn);
  csQn <- process(qn, verbose=log);
  print(csQn);
  stopifnot(identical(getNames(csQn), getNames(cs)));
  csQnList[[chipType]] <- csQn;
}





