library("aroma.affymetrix");
log <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the tiling array data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02", tags="Harvard,ROIs");
print(cdf);

csR <- AffymetrixCelSet$byName("MNtest", cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR, numChunks=20);
csM <- process(mn, verbose=more(log, 30));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get the "unique" CDF, which is generated if missing
cdfU <- getUniqueCdf(cdf, verbose=more(log, 60));
print(cdfU);

csU <- convertToUnique(csM, verbose=log);
print(csU);


dm <- matrix(1, ncol=2);
colnames(dm) <- sprintf("%s-%d", getNames(csU), 1:ncol(dm));
ms <- MatSmoothing(csU, design=dm);
print(ms);

dsMS <- process(ms, verbose=verbose);
print(dsMS);
