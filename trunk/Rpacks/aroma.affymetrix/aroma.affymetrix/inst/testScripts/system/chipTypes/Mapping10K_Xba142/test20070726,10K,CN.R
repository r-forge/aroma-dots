library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# Expected sample names
sampleNames <- c("0001-7", "0002-10", "0004-13", "0005-14", "0007-18", 
                      "0008-19", "0010-22", "2-DPrrr", "MH12", "MH18");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
cs <- extract(cs, keep);
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csAcc <- process(acc, verbose=log);
print(csAcc);
stopifnot(identical(getNames(csAcc), getNames(cs)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csAcc, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) Single array
theta <- extractMatrix(getFile(ces,1), verbose=log);
print(summary(theta));

data <- extractDataFrame(getFile(ces,1), addNames=TRUE, verbose=log);
print(summary(theta));

theta <- extractTheta(getFile(ces,1), drop=TRUE, verbose=log);
print(summary(theta));


# (b) Multiple arrays
theta <- extractMatrix(ces, verbose=log);
print(summary(theta));

data <- extractDataFrame(ces, addNames=TRUE, verbose=log);
print(summary(data));

theta <- extractTheta(ces, drop=TRUE, verbose=log);
print(summary(theta));




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);
cesFln <- process(fln, verbose=verbose);
print(cesFln);
stopifnot(identical(getNames(cesFln), getNames(ces)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Glad model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
glad <- GladModel(cesFln);
print(glad);
fit(glad, arrays=1, chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(glad);
print(ce);
process(ce, chromosomes=c(19,23), verbose=log);
## process(ce, verbose=log);
