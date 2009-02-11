###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson (and Henrik Bengtsson)
# Created: 2008-12-09
# Last modified: 2008-12-09
#
# Data set:
#  rawData/
#   MNtest/
#     Hs_PromPR_v02/
#       TESTIP.CEL
#       TESTIP.bar.txt (estimates by MAT)
###########################################################################
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
csM <- process(mn, verbose=more(log, 3));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=log);
print(csU);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to external estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare common units with prefix "chr1F".
cdf <- getCdf(csU);
units <- indexOf(cdf, "^chr1F");
cells <- getCellIndices(cdf, units=units, stratifyBy="pm", 
                                            unlist=TRUE, useNames=FALSE); 

# Get the chromosomal positions of these cells
acp <- AromaCellPositionFile$byChipType(getChipType(cdf));
pos <- acp[cells,2,drop=TRUE]; 

# Order cells by chromsomal position
o <- order(pos);
pos <- pos[o];
cells <- cells[o];

# Extract the corresponding signals of the first array
cf <- getFile(csU, 1);
y <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=log);


# Load external estimates
bar <- TabularTextFile("TESTIP.bar.txt", path=getPath(csR), 
                   columnNames=c("chromosome", "position", "signal"));
data <- readDataFrame(bar, colClasses=colClasses("cin"), nrow=435000);

# Extract the subset available in the aroma.affymetrix estimate
data <- subset(data, chromosome == "chr1" & position %in% pos);

# Order by position
o <- order(data$position);
data <- data[o,,drop=FALSE];

# Extract the external estimates
yB <- data$signal;

# Compare on the log2 scale
y <- log2(y);

stopifnot(length(yB) == length(y));

avgDiff <- mean(y-yB)^2;
cat(sprintf("avgDiff: %f\n", avgDiff));

plot(yB, y, pch=".");

stopifnot(avgDiff < 0.001);
