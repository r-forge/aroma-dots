library(aroma.affymetrix)

log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "HapMap270,100K,CEU,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get probe-sequence annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
apsList <- list();
for (chipType in chipTypes) {
  apsList[[chipType]] <- AromaProbeSequenceTextFile$byChipType(chipType);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRawList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
  stopifnot(identical(getNames(cs), sampleNames));
  csRawList[[chipType]] <- cs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csRawList;
csAccList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  acc <- AllelicCrosstalkCalibration(cs);
  print(acc);
  csC <- process(acc, verbose=log);
  print(csC);
  stopifnot(identical(getNames(csC), getNames(cs)));
  csAccList[[chipType]] <- csC;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-count normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csAccList;
csBcnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  bcn <- BaseCountNormalization(cs);
#  bcn <- BaseCountNormalization(cs, tags="*,lm", model="lm");
  print(bcn);
  csN <- process(bcn, verbose=log);
  print(csN);
  stopifnot(identical(getNames(csN), getNames(cs)));
  csBcnList[[chipType]] <- csN;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot probe sequence base count effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ff <- 3;
for (chipType in chipTypes) {
  aps <- apsList[[chipType]];
  counts <- countBases(aps, verbose=log);

  whats <- c("raw", "acc", "bcn");
  Ms <- list();
  for (what in whats) {
    if (what == "raw") {
      cs <- csRawList[[chipType]];
    } else if (what == "acc") {
      cs <- csAccList[[chipType]];
    } else if (what == "bcn") {
      cs <- csBcnList[[chipType]];
    }
    cf <- getFile(cs, ff);
    cfR <- getAverageFile(cs, verbose=log);
    yR <- readRawData(cfR, fields="intensities", drop=TRUE);
    y <- readRawData(cf, fields="intensities", drop=TRUE);
    M <- log2(y/yR);
    Ms[[what]] <- M;
    rm(y, yR, cfR, cs, M);
    gc <- gc();
  } # for (what ...)


  # Estimate standard deviations for log-ratios  
  print(sapply(Ms, FUN=function(M) { mad(diff(M), na.rm=TRUE) }));

  Mlim <- c(-1,1);
  Mlab <- expression(M == log[2](y/y[R]));
  xlim <- c(0, 25);
  for (what in names(Ms)) {
    M <- Ms[[what]];

    x11();
    layout(matrix(1:ncol(counts), nrow=2, byrow=TRUE));
    par(mar=c(5,4,2,1)+0.1);
    for (cc in 1:ncol(counts)) { 
      xlab <- sprintf("Number of %s:s", colnames(counts)[cc]);
      boxplot(M ~ counts[,cc], outline=FALSE, ylim=Mlim, xlim=xlim, ylab=Mlab, xlab=xlab);
      stext(side=3, pos=0, getFullName(cf));
      stext(side=3, pos=0.98, line=-1, cex=2, colnames(counts)[cc]);
      stextChipType(chipType);
    }
    rm(M);
  } # for (what ...)

  rm(counts);
  gc <- gc();
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csBcnList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- RmaCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
  print(plm);
  fit(plm, ram=1/2, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));
  cesCnList[[chipType]] <- ces;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesCnList <- cesCnList;
cesFlnList <- list();
for (chipType in names(csList)) {
  ces <- cesCnList[[chipType]];
  fln <- FragmentLengthNormalization(ces);
  print(fln);
  cesFln <- process(fln, verbose=verbose);
  print(cesFln);
  stopifnot(identical(getNames(cesFln), getNames(ces)));
  cesFlnList[[chipType]] <- cesFln;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Glad model test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Emulate list of ChipEffectSet:s where some arrays on exists in
# one of the sets
for (kk in seq(along=cesFlnList)) {
  ces <- cesFlnList[[kk]];
  ces <- extract(ces, setdiff(seq(ces), length(ces)+1-kk));
  cesFlnList[[kk]] <- ces;
}
glad <- GladModel(cesFlnList);
print(glad);

print(getTableOfArrays(glad));
nbrOfTestArrays <- nbrOfArrays(getSetTuple(glad));
nbrOfRefArrays <- nbrOfArrays(getReferenceSetTuple(glad));
stopifnot(identical(nbrOfTestArrays, nbrOfRefArrays));

fit(glad, arrays=1, chromosomes=19, verbose=log);

# Tests the case where one of the set does not have observations.
fit(glad, arrays=nbrOfArrays(glad), chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(glad);
print(ce);
process(ce, arrays=1:2, chromosomes=c(19,23), verbose=log);
