###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the SNPRMA 
# chip-effect estimates as estimated by oligo.
# It verifies that they give the same results whether or not one
# is normalizing towards the HapMap reference (as defined by oligo).
#
# Author: Henrik Bengtsson
# Created: 2008-12-04
# Last modified: 2008-12-05
#
# Data set:
#  rawData/
#   HapMap270,100K,CEU,testSet/
#     Mapping50K_Hind240/
#       NA06985,Hind,B5,3005533.CEL
#       NA06991,Hind,B6,3005533.CEL
#       NA06993,Hind,B4,4000092.CEL
#       NA06994,Hind,A7,3005533.CEL
#       NA07000,Hind,A8,3005533.CEL
#       NA07019,Hind,A12,4000092.CEL
###########################################################################

library("oligo");
library("aroma.affymetrix.wrappers");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

# Change this to FALSE if not normalizing toward the HapMap ref.
normalizeToHapmap <- TRUE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,6.0,CEU,testSet";
chipType <- "GenomeWideSNP_6";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet <- justSNPRMA(csR, normalizeToHapmap=normalizeToHapmap, verbose=log);
print(eSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet2 <- justSNPRMA(getPathnames(csR), normalizeToHapmap=normalizeToHapmap, verbose=TRUE);
print(eSet2);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare theta estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(eSet);
theta2 <- extractTheta(eSet2);

# Assert that the dimensions are the same
stopifnot(identical(dim(theta), dim(theta2)));

# Assert that the ordering of units and arrays are the same
stopifnot(identical(dimnames(theta), dimnames(theta2)));

# Assert that the estimates are very similar
tol <- 0.001;
stopifnot(all.equal(theta, theta2, tolerance=tol));
