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
plm <- RmaCnPlm(csAcc, mergeStrands=TRUE, combineAlleles=TRUE, 
                                              tags=c("+300", "*", "w"));
plm$shift <- +300;
plm$treatNAsAs <- "NA";
plm$treatNAsAs <- "weighted";
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
#  excludeChrXFromFit(fln);  # TO DO
print(fln);
cesFln <- process(fln, verbose=verbose);
print(cesFln);
stopifnot(identical(getNames(cesFln), getNames(ces)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GLAD model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
glad <- GladModel(cesFln);
print(glad);

fit(glad, arrays=1, chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesFln);
print(cbs);

fit(cbs, arrays=1, chromosomes=19, verbose=log);


csmList <- list(
  cbs  = CbsModel(cesFln),
  glad = GladModel(cesFln)
)

lapply(csmList, FUN=function(csm) {
  ce <- ChromosomeExplorer(csm);
  process(ce, arrays=1:3, chromosomes=c(1:2, 21:23), verbose=log);
})
