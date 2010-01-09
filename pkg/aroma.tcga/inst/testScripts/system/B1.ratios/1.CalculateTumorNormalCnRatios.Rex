############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Command line:
#  nohupR 1.PairedTumorNormalBAFs.Rex
#    --types=<types> --types=<centers> [--dataSets="<dataSet1>,...,<dataSetK>"]
#
# Description:
# Exports tumor-normal copy number ratios for samples of a given data set:
# 1. TCGA data sets are identified.
# 2. The user selects a data set.
# 3. For each tumor-normal pair, relative CNs are calculated and exported
#    to rawCnData/ adding a 'paired' tags to the data set.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Command line options
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!interactive()) {
  args <- commandArgs(asValues=TRUE, excludeReserved=TRUE, excludeEnvVars=TRUE)[-1];
  verbose && cat(verbose, "Raw command line arguments:");
  verbose && print(verbose, args);

  # Optional arguments
  if (!is.null(args$types)) {
    args$types <- trim(strsplit(args$types, split=", ")[[1]]);
  }
  if (!is.null(args$centers)) {
    args$centers <- trim(strsplit(args$centers, split=", ")[[1]]);
  }
  if (!is.null(args$dataSets)) {
    args$dataSets <- trim(strsplit(args$dataSets, split=", ")[[1]]);
  }
  if (!is.null(args$sampleNames)) {
    args$sampleNames <- trim(strsplit(args$sampleNames, split=", ")[[1]]);
  }
  verbose && cat(verbose, "Parsed command line arguments:");
  verbose && print(verbose, args);
}



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
} else {
  if (!is.null(args$types)) {
    uTypes <- intersect(uTypes, args$types);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  uCenters <- selectMenu(uCenters);
} else {
  if (!is.null(args$centers)) {
    uCenters <- intersect(uCenters, args$centers);
  }
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, uCenters);
dataSets <- dataSets[is.element(centers, uCenters)];

verbose && cat(verbose, "Filtered data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User selects a data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1)); 
} else {
  if (!is.null(args$dataSets)) {
    dataSets <- intersect(dataSets, args$dataSets);
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup tumor-normal paired allele B fraction data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsPairList <- lapply(dataSets, FUN=function(dataSet) {
  dsR <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", 
                                            paths="totalAndFracBData");
  setFullNamesTranslator(dsR, function(names, ...) {
    patterns <- BiospecimenCoreResource$getBarcodePatterns();
    pattern <- sprintf("^(%s)-(%s)[-]*(.*)", 
                       patterns$patient, patterns$sampleId);
    gsub(pattern, "\\1,\\5,\\8", names);
  });
  verbose && print(verbose, dsR);

  res <- extractTumorNormalPairs(dsR, verbose=verbose);
  verbose && print(verbose, res$tumor);
  verbose && print(verbose, res$normal);

  res;
});
verbose && print(verbose, dsPairList);

ns <- lapply(dsPairList, function(dsList) {
  length(dsList$tumor);
});
dsPairList <- dsPairList[ns > 0];
stopifnot(length(dsPairList) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Selecting samples to process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && cat(verbose, "Identified tumor-normal samples:");
sampleNames <- lapply(dsPairList, FUN=function(dsList) {
  getNames(dsList$tumor);
});
sampleNames <- sort(unique(unlist(sampleNames, use.names=FALSE)));
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);

# User can choose a subset of samples to process
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1),
                            title="Select samples of interest:");
  stopifnot(length(sampleNames) > 0);
} else {
  if (!is.null(args$sampleNames)) {
    sampleNames <- intersect(sampleNames, args$sampleNames);
  }
}

dsPairList <- lapply(dsPairList, FUN=function(dsList) {
  dsT <- extract(dsList$tumor, sampleNames, onMissing="drop");
  dsN <- extract(dsList$normal, sampleNames, onMissing="drop");
  # Sanity check
  stopifnot(identical(getNames(dsT), getNames(dsN)));
  list(tumor=dsT, normal=dsN);
});
ns <- lapply(dsPairList, function(dsList) {
  length(dsList$tumor);
});
dsPairList <- dsPairList[ns > 0];
stopifnot(length(dsPairList) > 0);

verbose && print(verbose, dsPairList);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate CN ratios, either as (i) M = log2(thetaT/thetaN), or 
# (ii) M = MT - MN, depending on the type of input data.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=dsPairList)) {
  dsPair <- dsPairList[[kk]];
  dsT <- dsPair$tumor;
  dsN <- dsPair$normal;
  stopifnot(getNames(dsT) == getNames(dsN));

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
} # for (kk ...)


############################################################################
# HISTORY:
# 2010-01-08
# o Added command line options.
# o Made the script to work on any number of data sets and chip types.
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
