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
# Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qn <- QuantileNormalization(cs);
print(qn);
csN <- process(qn, verbose=log);
print(csN);
stopifnot(identical(getNames(csN), getNames(cs)));


