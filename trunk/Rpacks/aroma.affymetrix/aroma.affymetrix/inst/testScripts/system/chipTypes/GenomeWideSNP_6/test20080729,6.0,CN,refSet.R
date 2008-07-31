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
cdf <- AffymetrixCdfFile$fromChipType(chipType, tags="Full");
print(cdf);

cs <- AffymetrixCelSet$fromName(dataSetName, cdf=cdf, verbose=log);
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


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

ces <- getChipEffectSet(plm);
print(ces);

fit(plm, verbose=log);


fln <- FragmentLengthNormalization(ces);
cesFln <- process(fln, verbose=log);
print(cesFln);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with specific reference set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use the robust average of the first three arrays as a reference
cesR <- extract(cesFln, 1:3);
ceR <- getAverageFile(cesR);
print(ceR);

sm <- CbsModel(cesFln, ceR);
print(sm);
fit(sm, arrays=1, chromosomes=19, verbose=log);
