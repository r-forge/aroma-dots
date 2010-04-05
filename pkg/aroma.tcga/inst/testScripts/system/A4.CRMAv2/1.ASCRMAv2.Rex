############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Command line:
#  nohupR 1.ASCRMAv2.Rex [--dataSets="<dataSet1>,...,<dataSetK>"]
#    [--sampleNames="<sampleName1>,...,<sampleNameK>"
#
# Description:
# 1. Identifies Affymetrix CEL data sets.
# 2. Runs ASCRMAv2 on these.
############################################################################
if (interactive()) savehistory();
library("aroma.affymetrix");
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

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

  # Optional arguments
  if (!is.null(args$dataSets)) {
    args$dataSets <- trim(strsplit(args$dataSets, split=", ")[[1]]);
  }
  if (!is.null(args$sampleNames)) {
    args$sampleNames <- trim(strsplit(args$sampleNames, split=", ")[[1]]);
  }
  verbose && cat(verbose, "Parsed command line arguments:");
  verbose && print(verbose, args);
}


rootPath <- "rawData/";
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all Affymetrix data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- "Genome_Wide_SNP_6";
dataSets <- grep(pattern, dataSets, value=TRUE);

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
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1));

  verbose && cat(verbose, "Filtered data sets:");
  verbose && print(verbose, dataSets);
} else {
  if (!is.null(args$dataSets)) {
    dataSets <- args$dataSets;
  }
}
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup list of data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipType <- "GenomeWideSNP_6,Full";
cdf <- AffymetrixCdfFile$byChipType(chipType);
verbose && print(verbose, cdf);

dsList <- list();
for (ii in seq(along=dataSets)) {
  dataSet <- dataSets[ii];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                            ii, dataSet, length(dataSets)));

  dsR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
  verbose && print(verbose, dsR);

  sdrfs <- SdrfFileSet$byPath(getPath(dsR));
  fnt <- TcgaFullNameTranslators$getBarcodeIdTranslator("sampleName,tumorType,rest");
  setFullNamesTranslator(dsR, list(sdrfs, fnt));

  verbose && print(verbose, dsR);
  # Nothing to do?
  if (nbrOfFiles(dsR) == 0) {
    verbose && cat(verbose, "No TCGA data files. Nothing to import.");
    verbose && exit(verbose);
    next;  
  }

  dsList[[ii]] <- dsR;
  
  verbose && exit(verbose);
} # for (ii ...)

# Drop empty
ns <- sapply(dsList, FUN=length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for all sample names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNamesList <- lapply(dsList, FUN=function(ds) {
  names <- getNames(ds);
  sort(unique(names));
});

sampleNames <- sort(unique(unlist(sampleNamesList, use.names=FALSE)));
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1));
} else {
  if (!is.null(args$sampleNames)) {
    sampleNames <- args$sampleNames;
  }
}
stopifnot(length(sampleNames) > 0);


# Drop files that do not contain the samples of interest
for (kk in seq(along=sampleNamesList)) {
  sampleNamesSet <- sampleNamesList[[kk]];
  ds <- dsList[[kk]];

  names <- getNames(ds);
  keep <- is.element(names, sampleNames);
  ds <- extract(ds, keep);

  dsList[[kk]] <- ds;
} # for (kk ...)

# Drop data sets that do not contain the samples of interest
ns <- sapply(dsList, FUN=length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);

verbose && cat(verbose, "Samples of interest:");
verbose && print(verbose, sampleNames);

## samplePatterns <- sprintf("^%s", sampleNames);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ASCRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ii in seq(along=dsList)) {
  dsR <- dsList[[ii]];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                      ii, getName(dsR), length(dsList)));
  verbose && print(verbose, dsR);

  # Import each data file
  verbose && enter(verbose, "ASCRMAv2");

  verbose && enter(verbose, "ASCRMAv2: Allelic crosstalk");
  acc <- AllelicCrosstalkCalibration(dsR, model="CRMAv2");
  verbose && print(verbose, acc);
  dsC <- process(acc, verbose=verbose)
  verbose && print(verbose, dsC);
  verbose && exit(verbose);

  verbose && enter(verbose, "ASCRMAv2: Nucleotide-position normalization");
  bpn <- BasePositionNormalization(dsC, target="zero");
  verbose && print(verbose, bpn);
  dsN <- process(bpn, verbose=verbose);
  verbose && print(verbose, dsN);
  verbose && exit(verbose);

  verbose && enter(verbose, "ASCRMAv2: Probe-level summarization");
  plm <- AvgCnPlm(dsN, mergeStrands=TRUE, combineAlleles=FALSE);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    fitCnProbes(plm, verbose=verbose);
    fit(plm, ram=5, verbose=verbose);
  }
  ces <- getChipEffectSet(plm)
  verbose && print(verbose, ces);
  verbose && exit(verbose);

  verbose && enter(verbose, "ASCRMAv2: Fragment-length normalization");
  fln <- FragmentLengthNormalization(ces, target="zero");
  verbose && print(verbose, fln);
  cesN <- process(fln, verbose=verbose);
  verbose && print(verbose, cesN);
  verbose && exit(verbose);
  
  verbose && exit(verbose);

  # Clean up
  rm(dsR, dsC, dsN, ces, cesN, acc, bpn, plm, fln);

  verbose && exit(verbose);
} # for (ii ...)



############################################################################
# HISTORY:
# 2010-04-05
# o Updated 'Usage'.
# 2010-01-18
# o Created.
############################################################################
