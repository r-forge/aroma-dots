############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Command line:
#  nohupR 1.MSCN,presmooting.Rex
#    --types=<types> --centers=<centers> --methods=<methods> [--dataSets="<dataSet1>,...,<dataSetK>"]
#
# Example:
#  nohupR testScripts/system/D2.MSCN/1.MSCN,presmooting.Rex --types=GBM --centers=broad.mit.edu --methods=ASCRMAv2,pairs &
#
# Description:
# 1. Finds all data sets with CN ratios.
# 2. Let user select by tumor type, centers, method, and data set names.
# 3. Runs the MSCN smoothing step on each remaining sample.
#
# Note that this script is optional, because MSCN will do it internally
# otherwise.  However, it is useful for distributed batch processing,
# because each data set can be smoothed independently of the others,
# e.g. as soon as they get produced.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("aroma.cn");
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
  if (!is.null(args$methods)) {
    args$methods <- trim(strsplit(args$methods, split=", ")[[1]]);
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
# Identify all tumor-normal CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPaths <- c("rawCnData", "cnData");
if (interactive()) {
  rootPath <- textMenu(rootPaths, value=TRUE);
} else {
  rootPath <- rootPaths[1];
}

rootPath <- Arguments$getReadablePath(rootPath, mustExist=TRUE);

patterns <- TcgaDccData$getDataSetPatterns();
pattern <- patterns$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,.*", patterns$dataset);
types <- gsub(pattern, "\\2", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select type:");
} else {
  if (!is.null(args$types)) {
    uTypes <- intersect(uTypes, args$types);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select preprocessing types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,", patterns$dataset);
types <- gsub(pattern, "", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select type:");
} else {
  if (!is.null(args$methods)) {
    uTypes <- intersect(uTypes, args$methods);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1)); 
} else {
  if (!is.null(args$dataSets)) {
    dataSets <- intersect(dataSets, args$dataSets);
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select samples to process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get all data sets
dsList <- lapply(dataSets, FUN=function(dataSet) {
  AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
});

sampleNames <- lapply(dsList, FUN=getNames);
sampleNames <- sort(unique(unlist(sampleNames, use.names=FALSE)));
stopifnot(length(sampleNames) > 0);

# User selection
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1));
} else {
  if (!is.null(args$sampleNames)) {
    sampleNames <- intersect(sampleNames, args$sampleNames);
  }
}
stopifnot(length(sampleNames) > 0);

# Keep only samples and data sets of interest
dsList <- lapply(dsList, extract, sampleNames, onMissing="drop");
ns <- sapply(dsList, length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);


dataSets <- unname(sapply(dsList, FUN=getFullName));
verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);


dataSets <- unname(sapply(dsList, FUN=getFullName));
verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

sampleNames <- lapply(dsList, FUN=getNames);
sampleNames <- sort(unique(unlist(sampleNames, use.names=FALSE)));
verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Target loci used for smoothing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb,HB20090503");
verbose && print(verbose, fitUgp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process sample by sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsSList <- list();
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                   kk, getFullName(ds), length(dsList)));

  sd <- 50e3;
  verbose && printf(verbose, "Kernel: %s\n", "Gaussian");
  verbose && printf(verbose, "Bandwidth (sd): %.2f\n", sd);

  sm <- TotalCnKernelSmoothing(ds, targetUgp=fitUgp, 
                                   kernel="gaussian", bandwidth=sd);
  verbose && print(verbose, sm);
  dsS <- process(sm, verbose=less(verbose, 1));
  verbose && print(verbose, dsS);

  dsSList[[kk]] <- dsS;

  verbose && exit(verbose);
} # for (kk ...)


############################################################################
# HISTORY:
# 2010-01-25
# o Now also MSCN normalized data can be smoothed.
# 2010-01-10
# o Created.
############################################################################
