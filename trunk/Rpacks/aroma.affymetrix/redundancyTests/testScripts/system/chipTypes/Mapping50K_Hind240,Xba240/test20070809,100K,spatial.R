library(aroma.affymetrix);
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

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
# Spatial intensity plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(csRawList);
setColorMaps(ae, "sqrt,yellow");
print(ae);
if (interactive())
  process(ae, verbose=log);
stopifnot(identical(unname(getArrays(ae)), getNames(cs)));
