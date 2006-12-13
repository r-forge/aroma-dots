library(oligo);
library(aroma.affymetrix);
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
.Machine$float.eps <- sqrt(.Machine$double.eps);
tolerance <- .Machine$float.eps^0.5;

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

tests <- c("pd", "raw", "pcr", "qn", "plm", "snpQSet");

#tests <- c("plm");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: Raw data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("pd" %in% tests) {
  validate(pd, subset=1/100, verbose=verbose);
}

if (!exists("snpFSet")) {
  snpFSet <- extractSnpFeatureSet(ds, verbose=verbose);
  verbose && print(verbose, snpFSet);
}

if ("raw" %in% tests) {
  verbose && enter(verbose, "Verifying reading of raw data");
  # Select random cells to validate
  cells <- sample(nbrOfCells(cdf), size=1000);

  # Read data using aroma.affymetrix
  verbose && enter(verbose, "Using aroma.affymetrix");
  dataA <- readCelIntensities(getPathnames(ds), indices=cells);
  dimnames(dataA) <- NULL;
  verbose && exit(verbose);
  
  # Read data using oligo
  verbose && enter(verbose, "Using oligo");
  dataO <- exprs(snpFSet)[map[cells],];
  dimnames(dataO) <- NULL;
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
  sc <- PcrProbeCalibration(ds, flavor="oligo");
  dsC <- process(sc, verbose=verbose);
  verbose && exit(verbose);

  # Calibrate using oligo
  verbose && enter(verbose, "Using oligo");
  if (!exists("pmsC")) {
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
    rm(pms);
    pmsC <- tmp[[1]];
    rm(SeqMat, L, tmp);
  }
  verbose && exit(verbose);

  # Comparing
  verbose && enter(verbose, "Comparing");
  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  dataA <- readCelIntensities(getPathnames(dsC), indices=cells);
  dimnames(dataA) <- NULL;

  rr <- match(map[cells], as.integer(rownames(pmsC)));
  dataO <- pmsC[rr,];
  dimnames(dataO) <- NULL;
  stopifnot(all.equal(dataA, dataO, tolerance=tolerance));
  verbose && exit(verbose);

  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("qn" %in% tests) {
  verbose && enter(verbose, "Verifying quantile normalization");

  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  # Normalize using aroma.affymetrix
  verbose && enter(verbose, "Using aroma.affymetrix");
  norm <- OligoQuantileNormalization(dsC);
  verbose && print(verbose, norm);
  dsN <- process(norm, verbose=verbose);
  verbose && print(verbose, dsN);
  verbose && exit(verbose);

  # Normalize using oligo
  verbose && enter(verbose, "Using oligo");
  if (!exists("pmsN")) {
    yTarget <- getReferenceQuantiles(pd);
    pmsN <- oligo::normalizeToSample(pmsC, yTarget);
    rm(yTarget);
    rm(snpFSetN);
  }
  if (!exists("snpFSetN")) {
    snpFSetN <- snpFSet;
    foo <- oligo::"pm<-";
    snpFSetN <- foo(snpFSetN, pmsN);
  }
  verbose && exit(verbose);

  # Comparing
  verbose && enter(verbose, "Comparing");
  # Select random SNPs to validate
  units <- sample(indexOf(cdf, pattern="^SNP"), size=1000);
  cells <- unlist(getCellIndices(cdf, units=units, stratifyBy="pm"), use.names=FALSE);

  dataA <- readCelIntensities(getPathnames(dsN), indices=cells);
  dimnames(dataA) <- NULL;

  rr <- match(map[cells], as.integer(rownames(pmsC)));
  dataO <- pmsN[rr,];
  dimnames(dataO) <- NULL;

  stopifnot(all.equal(dataA, dataO, tolerance=tolerance));
  verbose && exit(verbose);

  verbose && exit(verbose);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: Quantile normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("plm" %in% tests) {
  verbose && enter(verbose, "Verifying PLM");

  verbose && enter(verbose, "Using aroma.affymetrix");
  plm <- RmaSnpPlm(dsN, mergeStrands=FALSE, flavor="oligo");
  verbose && print(verbose, plm);
  unitsDone <- fit(plm, moreUnits=moreUnits, verbose=verbose);
  ces <- getChipEffects(plm);
  verbose && exit(verbose);

  
  verbose && enter(verbose, "Using oligo");
  if (!exists("theta")) {
    ltheta <- summSnp(snpFSetN, method=3, normalize=FALSE, background=FALSE);
    theta <- 2^ltheta;
    rm(snpQSet);
  }
  if (!exists("snpQSet")) {
    snpQSet <- new("SnpQSet",
            senseThetaA=ltheta["sense","A",,],
            senseThetaB=ltheta["sense","B",,],
            antisenseThetaA=ltheta["antisense","A",,],
            antisenseThetaB=ltheta["antisense","B",,],
            phenoData=phenoData(snpFSetN),
            experimentData=experimentData(snpFSetN),
            annotation=annotation(snpFSetN));
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Comparing 'raw' chip effects");
  units <- sample(indexOf(cdf, pattern="^SNP"), size=500);
  unitNames <- getUnitNames(cdf, units=units);
  gs <- readCdfGroupStrands(getPathname(cdf), units=units, what="target");
  cesU <- ces[units];
  thetaU <- theta[,,unitNames,];

  for (kk in seq(along=units)) {
    if (kk %% 100 == 0)
      cat(kk, ",", sep="");

    dataA <- cesU[[kk]];
    dataA <- list(
      matrix(unlist(dataA[1:2]), nrow=2, byrow=TRUE),
      matrix(unlist(dataA[3:4]), nrow=2, byrow=TRUE)
    )
    names(dataA) <- matrix(gs[[kk]], nrow=2)[1,];
    o <- order(names(dataA));
    dataA <- dataA[o];

    dataO <- thetaU[,,kk,,drop=TRUE];
    dataO <- list(antisense=dataO["antisense",,], sense=dataO["sense",,]);
    dataO <- lapply(dataO, FUN=function(X) { dimnames(X) <- NULL; X; })

    stopifnot(all.equal(dataA, dataO, tolerance=tolerance));
  }
  cat("\n");
  verbose && exit(verbose);

  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TEST: extractSnpQSet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ("snpQSet" %in% tests) {
  verbose && enter(verbose, "Verifying extractSnpQSet()");
  if (!exists("snpQSet2")) {
    snpQSet2 <- extractSnpQSet(ces, verbose=verbose);
  }
  
  for (fcn in list(senseThetaA, senseThetaB, 
                                 antisenseThetaA, antisenseThetaB)) {
    stopifnot(all.equal(fcn(snpQSet), fcn(snpQSet2), tolerance=tolerance));
  }
  verbose && exit(verbose);
}
