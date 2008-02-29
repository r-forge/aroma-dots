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
cdf <- AffymetrixCdfFile$fromChipType(chipType, tags="Full,r2");
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
process(ae, verbose=log);
