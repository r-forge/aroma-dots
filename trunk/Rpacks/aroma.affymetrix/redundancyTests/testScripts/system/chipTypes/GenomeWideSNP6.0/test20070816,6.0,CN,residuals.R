library(aroma.affymetrix)

log <- Arguments$getVerbose(-50);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap270,100K,CEU,testSet";
chipType <- "GenomeWideSNP_6";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(csC);
print(plm);
fit(plm, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Residuals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rs <- calculateResidualSet(plm, verbose=log);
print(rs);

ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);
if (interactive())
  process(ae, verbose=log);
