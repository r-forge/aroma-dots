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
# Fitting log-additive model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(cs);
print(plm);
fit(plm, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Some basic quality scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- getChipEffectSet(plm);

# Boxplots of via log2(theta), RLE, and NUSE
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplot(ces, type="theta", transform=log2);
plotBoxplot(ces, type="RLE", arrays=c(2,4:6));
plotBoxplot(ces, type="NUSE");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Advanced usage
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculating statistics without plotting
theta <- boxplotStats(ces, type="theta", transform=log2);
nuse <- boxplotStats(ces, type="NUSE");
rle <- boxplotStats(ces, type="RLE");

# Subset of arrays (avoids calculating stats for all arrays)
rleB <- boxplotStats(ces, type="RLE", arrays=c(2,4:6));
for (name in names(rleB)) {
  stopifnot(identical(rleB[name], rle[name]));
}

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotBoxplotStats(theta, main="theta");
plotBoxplotStats(rle[c(2,4:6)], main="RLE");
plotBoxplotStats(nuse, main="NUSE");


# Calculates unit-specific RLE and NUSE scores
units <- 1000+1:5000;
theta <- extractMatrix(ces, units=units);
rle <- extractMatrix(ces, units=units, field="RLE");
nuse <- extractMatrix(ces, units=units, field="NUSE");

# Plotting the above statistics
layout(matrix(1:4, ncol=2, byrow=TRUE));
plotDensity(log2(theta), main="theta");
plotDensity(rle[,c(2,4:6)], main="RLE");
plotDensity(nuse, main="NUSE");

# ...same, but basic unit annotation data added
units <- 1000+1:500;
theta <- extractDataFrame(ces, units=units, addNames=TRUE);
rle <- extractDataFrame(ces, units=units, field="RLE", addNames=TRUE);
nuse <- extractDataFrame(ces, units=units, field="NUSE", addNames=TRUE);
