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
  csAcc <- process(acc, verbose=log);
  print(csAcc);
  stopifnot(identical(getNames(csAcc), getNames(cs)));
  csAccList[[chipType]] <- csAcc;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csAccList;
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
# Robust estimate of standard deviation of raw CNs
# (The default is to use a first-order difference variance estimator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate it for some autosomal chromosomes and ChrX
cns <- CbsModel(cesFlnList);
res <- estimateSds(cns, chromosomes=c(1:6, 23), verbose=log);
chrX <- which(rownames(res) == "23");

layout(matrix(1:2, nrow=2));
par(mar=c(4,4,0.5,2)+0.1);

xlim <- c(1,(ncol(res)+2));
ylim <- c(0, 1.05*max(res, na.rm=TRUE));
ylab <- expression(hat(sigma)==s[Delta](log2(theta/theta[R])));
cols <- seq(length=nrow(res));
ltys <- rep(1, times=nrow(res));
cols[chrX] <- "blue";
ltys[chrX] <- 4;

plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab=ylab);
for (kk in seq(length=nrow(res))) {
  chromosome <- rownames(res)[kk];
  points(res[kk,], col=cols[kk], pch=19);
  lines(res[kk,], col=cols[kk], lty=ltys[kk], lwd=2);
}
legend("topright", col=cols, lty=ltys, pch=19, lwd=2, 
                   legend=sprintf("Chr%s", rownames(res)), bty="n");


df <- as.data.frame(res[-chrX,,drop=FALSE]);
colnames(df) <- seq(length=nrow(df));
boxplot(df, ylim=ylim, ylab=ylab, xlab="Array");
points(res[chrX,], col="blue", pch=19, cex=1.5);
legend("bottomright", col=c("black", cols[chrX]), pch=19, lwd=2, 
                   legend=c("Autosomal", "ChrX"), horiz=TRUE, bty="n");
devDone();
