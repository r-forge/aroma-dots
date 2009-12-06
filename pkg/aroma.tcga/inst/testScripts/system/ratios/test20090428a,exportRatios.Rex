############################################################################
# Author: Henrik Bengtsson, hb@stat.berkeley.edu
#
# Description:
# Exports tumor-normal copy number ratios for samples of a given data set.
# All TCGA ovarian data sets are identified.  The user choose which to
# export.
############################################################################
if (interactive()) savehistory();
library("aroma.core");
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User selects: data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

dataSet <- textMenu(dataSets, value=TRUE); 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup (total,fracB) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- NULL;
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*", paths="totalAndFracBData");
setFullNamesTranslator(ds, function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
});
verbose && print(verbose, ds);

# Extract tuples of arrays for each unique sample
sampleNames <- unique(getNames(ds));
idxs <- indexOf(ds, sampleNames);
if (!is.matrix(idxs)) {
  # Drop non-pairs
  ns <- sapply(idxs, length);
  verbose && print(verbose, table(ns));
  isPair <- (ns == 2);
  idxs <- indexOf(ds, sampleNames[isPair]);
  stopifnot(is.matrix(idxs));
}
stopifnot(nrow(idxs) == 2);
dimnames(idxs) <- names(idxs) <- NULL;
colnames(idxs) <- getNames(ds)[idxs[1,]];

ds2 <- extract(ds, idxs);
idxs <- matrix(seq(ds2), nrow=2);

# Identify sample type, e.g. tumor or normal etc.
types <- sapply(ds2, function(df) getTags(df)[1]);
types <- gsub("[A-Z]$", "", types);
dim(types) <- dim(idxs);
uTypes <- unique(as.vector(types));
verbose && cat(verbose, "Identified types:");
verbose && print(verbose, uTypes);

# Sort each column according to sample type
o <- apply(types, MARGIN=2, FUN=order);
for (kk in seq(length=ncol(o))) {
  types[,kk] <- types[o[,kk],kk];
  idxs[,kk] <- idxs[o[,kk],kk];
}
stopifnot(all(is.finite(idxs)));

dsT <- extract(ds2, idxs[1,]);
dsN <- extract(ds2, idxs[2,]);
verbose && cat(verbose, "Data set of tumors:");
verbose && print(verbose, dsT);
verbose && cat(verbose, "Data set of normals:");
verbose && print(verbose, dsN);

# Sanity check
stopifnot(identical(getNames(dsT), getNames(dsN)));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate CN ratios, either as (i) M = log2(thetaT/thetaN), or 
# (ii) M = MT - MN, depending on the type of input data.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsC <- exportTotalCnRatioSet(dsT, ref=dsN, logBase=2, tags="pairs", 
                                                      verbose=verbose);
setFullNamesTranslator(dsC, function(names, ...) {
  barcodePattern <- "(TCGA-[0-9]{2}-[0-9]{4}),([0-9]{2}[A-Z]),([^,]*)";
  pattern <- sprintf("^%s,ref=%s,(.*)$", barcodePattern, barcodePattern);
  gsub(pattern, "\\1,\\2vs\\5,\\3vs\\6,\\7", names);
});
verbose && cat(verbose, "Data set of tumor-normal copy numbers:");
verbose && print(verbose, dsC);
verbose && cat(verbose, "All full names:");
verbose && print(verbose, getFullNames(dsC));

stop();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsC <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=c(tags, "pairs"), chipType="*");
setFullNamesTranslator(dsC, function(names, ...) {
  barcodePattern <- "(TCGA-[0-9]{2}-[0-9]{4}),([0-9]{2}[A-Z]),([^,]*)";
  pattern <- sprintf("^%s,ref=%s,(.*)$", barcodePattern, barcodePattern);
  gsub(pattern, "\\1,\\2vs\\5,\\3vs\\6,\\7", names);
});
print(dsC);
print(getFullNames(dsC));

aa <- indexOf(dsC, "TCGA-09-0369");
df <- getFile(dsC, aa);
name <- paste(c(getName(df), getTags(df)[1]), collapse=",");
cn <- extractRawCopyNumbers(df, chromosome=1);
plot(cn, pch=".", xScale=1e-6);
stext(side=3, pos=0, name);


############################################################################
# HISTORY:
# 2009-09-24
# o Added header comments.
# o Added verbose output.
# 2009-05-04
# o Now utilizing R.menu for selecting data set.
# 2009-04-24
# o Created.
############################################################################
