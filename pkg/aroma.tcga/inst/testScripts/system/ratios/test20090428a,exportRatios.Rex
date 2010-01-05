############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/ratios/");
#
# Description:
# Exports tumor-normal copy number ratios for samples of a given data set:
# 1. TCGA data sets are identified.
# 2. The user selects a data set.
# 3. For each tumor-normal pair, relative CNs are calculated and exported
#    to rawCnData/ adding a 'paired' tags to the data set.
############################################################################
if (interactive()) savehistory();
library("aroma.core");
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getTypeTag <- function(this, ...) {
  gsub("[A-Z]$", "", getTags(df)[1]);
}


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
ds0 <- ds;

# Extract tuples of arrays for each unique sample
sampleNames <- unique(getNames(ds));
idxs <- indexOf(ds, sampleNames);
idxs <- split(idxs, names(idxs));
ns <- sapply(idxs, length);
verbose && print(verbose, table(ns));

# Drop tuples for which lack either a tumor or a normal
types <- sapply(ds, getTypeTag);
idxs <- lapply(idxs, FUN=function(idxs) {
  types <- types[idxs];
  tIdx <- which(is.element(types, c("01"      )))[1];
  nIdx <- which(is.element(types, c("10", "11")))[1];
  idxs <- idxs[c(tIdx, nIdx)];
  idxs;  
});

keep <- sapply(idxs, FUN=function(x) all(is.finite(x)));
idxs <- idxs[keep];

# Sanity check
ns <- sapply(idxs, length);
stopifnot(all(ns == 2));

dim <- c(2, length(idxs));

ds <- extract(ds, unlist(idxs));

# Sanity checks
types <- sapply(ds, getTypeTag);
dim(types) <- dim;
stopifnot(all(types[1,] == "01"));
stopifnot(all(is.element(types[2,], c("10", "11"))));

# Sanity checks
names <- sapply(ds, getName);
stopifnot(all(names[1,] == names[2,]));

# Create an index matrix
idxs <- seq(ds);
dim(idxs) <- dim;
rownames(idxs) <- c("T", "N");
colnames(idxs) <- names[1,];

uTypes <- unique(as.vector(types));
verbose && cat(verbose, "Identified types:");
verbose && print(verbose, uTypes);

# Extract tumors and normals
dsT <- extract(ds, idxs["T",]);
dsN <- extract(ds, idxs["N",]);
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


############################################################################
# HISTORY:
# 2010-01-04
# o Before the script would drop any samples with (n != 2) replicates.
#   Now the script picks out a tumor-normal pair also for n > 2.
# o Updated to work with new R.filesets v0.7.0, whose indexOf() no longer
#   returns a list if there are multiple matches.
# 2009-09-24
# o Added header comments.
# o Added verbose output.
# 2009-05-04
# o Now utilizing R.menu for selecting data set.
# 2009-04-24
# o Created.
############################################################################
