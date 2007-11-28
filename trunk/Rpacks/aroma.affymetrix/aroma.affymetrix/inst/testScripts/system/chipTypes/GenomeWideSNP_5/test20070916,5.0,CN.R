library(aroma.affymetrix)

log <- Verbose(threshold=-50, timestamp=TRUE);

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
# Identify cells *not* for chromosome X & Y
cdf <- getCdf(cs);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosome(gi, 23:24);
cells <- getCellIndices(cdf, units=units, useNames=FALSE, unlist=TRUE);
rm(units);
cells <- setdiff(1:nbrOfCells(cdf), cells);
acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=cells, tags=c("*", "-XY"));
rm(cells);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Averaging probe-level model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE);
print(plm);
fit(plm, verbose=log);


ces <- getChipEffectSet(plm);
theta <- extractMatrix(ces, units=1000:1002);

fln <- FragmentLengthNormalization(ces);
#cesFln <- process(fln, verbose=log);
cesFln <- ces;
print(cesFln);

cnr <- CbsModel(cesFln);
print(cnr);

ce <- ChromosomeExplorer(cnr);
print(ce);
process(ce, arrays=1, chromosomes=19:23, verbose=log);
