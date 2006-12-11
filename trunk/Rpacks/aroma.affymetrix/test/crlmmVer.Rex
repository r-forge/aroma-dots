library(oligo);
library(aroma.affymetrix);
source("init.R");

.Machine$float.eps <- sqrt(.Machine$double.eps);

moreUnits <- 30;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "raw/Test500K/Mapping250K_Nsp";

path <- "raw/Test100K/Mapping50K_Hind240";

ds <- AffymetrixCelSet$fromFiles(path);
print(ds);

cdf <- getCdf(ds);
pd <- PlatformDesign(cdf);

# CDF-to-PD index map
map <- getReadMap(pd);

tests <- c("pd", "raw", "pcr", "qn");

#tests <- c("qn");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: Raw data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("pd" %in% tests) {
  validate(pd, subset=1/100, verbose=verbose);
}

if (!exists("snpFSet")) {
  snpFSet <- extractSnpFeatureSet(ds, verbose=verbose);
}

if ("raw" %in% tests) {
  verbose && enter(verbose, "Verifying reading of raw data");
  # Select random cells to validate
  cells <- sample(nbrOfCells(cdf), size=1000);

  # Read data using aroma.affymetrix
  verbose && enter(verbose, "Using aroma.affymetrix");
  dataA <- readCelIntensities(getPathnames(ds), indices=cells);
  colnames(dataA) <- basename(colnames(dataA));
  verbose && exit(verbose);
  
  # Read data using oligo
  verbose && enter(verbose, "Using oligo");
  dataO <- exprs(snpFSet)[map[cells],];
  rownames(dataO) <- NULL;
  verbose && exit(verbose);

  stopifnot(identical(dataA, dataO));
  verbose && exit(verbose);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: PCR calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("pcr" %in% tests) {
  verbose && enter(verbose, "Verifying calibration for PCR effects");

  # Calibrate using aroma.affymetrix
  verbose && enter(verbose, "Using aroma.affymetrix");
  sc <- PcrProbeCalibration(ds, subsetToFit="asOligo");
  dsC <- process(sc, verbose=verbose);
  verbose && exit(verbose);

  # Calibrate using oligo
  verbose && enter(verbose, "Using oligo");
  pns <- oligo::probeNames(snpFSet);
  pms <- oligo::pm(snpFSet);
  verbose && enter(verbose, "Creating design matrix");
  SeqMat <- sequenceDesignMatrix(pmSequence(snpFSet));
  data(list=annotation(snpFSet));
  theLengths <- annot[match(pns, annot$SNP), "Length"];
  med <- median(theLengths, na.rm=TRUE);
  theLengths[is.na(theLengths)] <- med;
  L <- ns(theLengths, df=3);
  rm(med, theLengths);
  verbose && exit(verbose);
  verbose && enter(verbose, "Calibrating");
  tmp <- CorrectSequenceLength(pms, cbind(SeqMat, L));
  pms <- tmp[[1]];
  rm(SeqMat, L, tmp);
  verbose && exit(verbose);

  # Comparing
  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  dataA <- readCelIntensities(getPathnames(dsC), indices=cells);
  dimnames(dataA) <- NULL;

  rr <- match(map[cells], as.integer(rownames(pms)));
  dataO <- pms[rr,];
  dimnames(dataO) <- NULL;
  stopifnot(all.equal(dataA, dataO, tolerance=.Machine$float.eps^0.5));
  verbose && exit(verbose);

  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("qn" %in% tests) {
  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  # Normalize using aroma.affymetrix
  verbose && enter(verbose, "Using aroma.affymetrix");
  hapmapQs <- getReferenceQuantiles(pd);
#  norm <- OligoQuantileNormalization(dsC);
  norm <- QuantileNormalization(dsC, typesToUpdate="pm", 
               targetDistribution=hapmapQs, tags=c("*", "HapMap"));
  verbose && print(verbose, norm);
  dsN <- process(norm, verbose=verbose);
  verbose && print(verbose, dsN);
  verbose && exit(verbose);

  # Normalize using oligo
  verbose && enter(verbose, "Using oligo");
  data(list=paste(platform(snpFSet), "Ref", sep=""));
  pmsN <- normalizeToSample(pms, reference);
  rm(reference);
  verbose && exit(verbose);

  # Comparing
  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  dataA <- readCelIntensities(getPathnames(dsN), indices=cells);
  dimnames(dataA) <- NULL;

  rr <- match(map[cells], as.integer(rownames(pms)));
  dataO <- pms[rr,];
  dimnames(dataO) <- NULL;

  # Check correlation
  cor <- diag(cor(dataA, dataO));
  verbose, print(verbose, cor);
  if (any(cor < 0.95)) {
    throw("Correlation between aroma.affymetrix and oligo is < 0.95");
  }
}

stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quantile normalization of PMs toward HapMap references
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit RMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaSnpPlm(dsN, mergeStrands=FALSE);
verbose && print(verbose, plm);
units <- fit(plm, moreUnits=moreUnits, verbose=verbose);
verbose && exit(verbose);

