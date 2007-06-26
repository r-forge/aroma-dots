library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "Affymetrix_2004-100k_trios,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a CEL set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipTypes[1], verbose=log);
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set a custom CDF (here just use the other CDF of the chip in the set)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$fromChipType(chipTypes[2]);
setCdf(cs, cdf);
print(cs);
stopifnot(identical(getCdf(cs), cdf));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA background correction with custom CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bg <- RmaBackgroundCorrection(cs);
print(bg);
csBg <- process(bg, verbose=log);
print(csBg);
stopifnot(identical(getCdf(csBg), cdf));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quantile normalization with custom CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qn <- QuantileNormalization(cs);
print(qn);
csQn <- process(qn, verbose=log);
print(csQn);
stopifnot(identical(getCdf(csQn), cdf));

