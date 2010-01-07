############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
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
# User selects: data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patterns <- TcgaDccData$getDataSetPatterns();
pattern <- patterns$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s.*", patterns$dataset);
types <- gsub(pattern, "\\2", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1),
                       title="Select type:");
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User selects a data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSet <- textMenu(dataSets, value=TRUE); 
} else {
  dataSet <- dataSets[1];
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup (total,fracB) data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- NULL;
dsR <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*", paths="totalAndFracBData");
setFullNamesTranslator(dsR, function(names, ...) {
  patterns <- BiospecimenCoreResource$getBarcodePatterns();
  pattern <- sprintf("^(%s)-(%s)[-]*(.*)", 
                     patterns$patient, patterns$sampleId);
  gsub(pattern, "\\1,\\5,\\8", names);
});
verbose && print(verbose, dsR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract all tumors
dsT <- extractByTcgaType(dsR, pattern="^01");

# Extract all normals
dsN <- extractByTcgaType(dsR, pattern="^1[01]");

# Keep only samples for which there is both a tumor and a normal
sampleNames <- intersect(getNames(dsT), getNames(dsN));
sampleNames <- unique(sampleNames);
verbose && cat(verbose, "Identified tumor-normal samples:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);

# User can choose a subset of samples to process
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1),
                            title="Select samples to normalize:");
}
stopifnot(length(sampleNames) > 0);

dsT <- extract(dsT, sampleNames);
dsN <- extract(dsN, sampleNames);

# Drop duplicated samples
dsT <- extract(dsT, !duplicated(getNames(dsT)));
dsN <- extract(dsN, !duplicated(getNames(dsN)));

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
