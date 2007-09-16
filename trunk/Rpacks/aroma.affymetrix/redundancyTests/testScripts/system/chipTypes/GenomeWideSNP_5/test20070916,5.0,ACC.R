library(aroma.affymetrix)

log <- Arguments$getVerbose(-50);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap270,5.0,CEU,testSet";
chipType <- "GenomeWideSNP_5";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993",
                 "NA07019", "NA07022", "NA07056");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL set and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
print(cs);
stopifnot(identical(getNames(cs), sampleNames));
cdf <- getCdf(cs);
print(cdf);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify cells *not* for chromosome X
cdf <- getCdf(cs);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosome(gi, 23);
cells <- getCellIndices(cdf, units=units, useNames=FALSE, unlist=TRUE);
rm(units);
cells <- setdiff(1:nbrOfCells(cdf), cells);
acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=cells, tags=c("*", "-X"));
rm(cells);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));
