library(aroma.affymetrix)

log <- Verbose(threshold=-50, timestamp=TRUE);

dataSetName <- "HapMap270,6.0,CEU,testSet";
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
acc <- AllelicCrosstalkCalibration(cs);
print(acc);

csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Averaging probe-level model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);
fit(plm, verbose=log);


ces <- getChipEffectSet(plm);
theta <- extractMatrix(ces, units=1000:1002);

fln <- FragmentLengthNormalization(ces);
cesFln <- process(fln, verbose=log);
print(cesFln);

cnr <- CbsModel(cesFln);
print(cnr);

ce <- ChromosomeExplorer(cnr);
print(ce);
process(ce, arrays=1, chromosomes=c(19,23), verbose=log);
