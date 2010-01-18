############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/normalize/");
#
# Command line:
#  nohupR 1.PrincipalCurveNormalization.Rex
#    --types=<types> --centers=<centers> --methods=<methods> 
#    [--dataSets="<dataSet1>,...,<dataSetK>"]
#
# Description:
# Normalized (tumor,normal) CN ratios using principal-curve normalization.
#  1. ...
#  2. ...
############################################################################
if (interactive()) savehistory();
library("aroma.cn");
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tcgaFNT01 <- function(names, ...) {
  patterns <- BiospecimenCoreResource$getBarcodePatterns();
  pattern <- sprintf("^(%s)-(%s)[-]*(.*)", 
                     patterns$patient, patterns$sampleId);
  gsub(pattern, "\\1,\\5,\\8", names);
}


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

# Drop already PCN normalized ones
keep <- (regexpr(",PCN", dataSets) == -1);
dataSets <- dataSets[keep];

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
verbose && print(verbose, dataSets);
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select preprocessing types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,", patterns$dataset);
types <- gsub(pattern, "", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select preprocessing methods:");
} else {
  if (!is.null(args$methods)) {
    uTypes <- intersect(uTypes, args$methods);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);

# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1), title="Select dataset(s):");
}
stopifnot(length(dataSets) > 0);


verbose && cat(verbose, "Selected data sets:");
verbose && print(verbose, dataSets);


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
}

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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit on log scale or not?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  choises <- c(TRUE, FALSE);
  onLog2Scale <- textMenu(choises, title="Fit on log2 scale?", value=TRUE); 
} else {
  onLog2Scale <- FALSE;
}
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process data set by data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=dsList)) {
  ds <- dsList[[kk]];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                 kk, getFullName(ds), length(dsList)));

  setFullNamesTranslator(ds, tcgaFNT01);
  
  verbose && cat(verbose, "Data set:");
  verbose && print(verbose, ds);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract tumor-normal pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract all tumors
  dsT <- extractByTcgaType(ds, pattern="^(01|20)");
  
  # Extract all normals
  dsN <- extractByTcgaType(ds, pattern="^(10|11)");
  
  # Keep only samples for which there is both a tumor and a normal
  sampleNames <- intersect(getNames(dsT), getNames(dsN));
  sampleNames <- unique(sampleNames);
  verbose && cat(verbose, "Identified tumor-normal samples:");
  verbose && print(verbose, sampleNames);

  # Nothing to do?
  if (length(sampleNames) == 0) {
    verbose && cat(verbose, "No paired tumor-normals. Skipping.");
    verbose && exit(verbose);
  }
  
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
  # Principal curve normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setFullNamesTranslator(dsT, NULL);
  setFullNamesTranslator(dsN, NULL);
  pcn <- PrincipalCurveNormalization(dsT, dsN);
  # Fit on the log2 scale?
  if (onLog2Scale) {
    pcn$.transformFcn <- function(x) log2(x);
    pcn$.untransformFcn <- function(y) 2^y;
    setTags(pcn, c("*", "log2"));
  }
  verbose && print(verbose, pcn);
  
  dsTN <- process(pcn, verbose=verbose);
  
  verbose && cat(verbose, "Normalized data set:");
  verbose && print(verbose, dsTN);

  verbose && exit(verbose);
} # for (kk ...)


############################################################################
# HISTORY:
# 2010-01-17
# o Extended to process multiple data sets
# 2010-01-05
# o Updated the plot cut'n'paste code at the end.
# o Updated script so that one can fit on the log2-scale (default).
# o Updated code for identifying tumor-normal pairs.
# o Added instructions to the header.
# 2009-07-15
# o Created.
############################################################################
