library(aroma.affymetrix)

log <- Arguments$getVerbose(-50);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap6.0-testSet";
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
# Spatial intensity plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(cs);
setColorMaps(ae, c("sqrt,yellow", "log2,yellow"));
print(ae);
process(ae, verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qn <- QuantileNormalization(cs);
print(qn);
csN <- process(qn, verbose=log);
print(csN);
stopifnot(identical(getNames(csN), getNames(cs)));


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
plm <- RmaPlm(csN);
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
