############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/import/");
#
# Description:
# 1. Identifies CRMAv2 (aroma.affymetrix) chip-effect data sets.
# 2. Let user choose samples to export by tumor type and data set(s).
# 3. Exports selected *,chipEffects.CEL data files to totalAndFracBData/
#    where the data set tag is 'ASCRMAv2'.
############################################################################
if (interactive()) savehistory();
library("aroma.affymetrix");
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

rootPath <- "plmData/";
rootPath <- Arguments$getReadablePath(rootPath);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);

verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all CRMA v2 data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
crmaTags <- c(ASCRMAv2="ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
pattern <- sprintf(".*,%s$", crmaTags);
dataSets <- grep(pattern, dataSets, value=TRUE);

verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  center <- textMenu(uCenters, value=TRUE);
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, center);
dataSets <- dataSets[centers == center];

verbose && cat(verbose, "Filtered data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1));

  verbose && cat(verbose, "Filtered data sets:");
  verbose && print(verbose, dataSets);
  stopifnot(length(dataSets) > 0);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup list of data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list();
for (ii in seq(along=dataSets)) {
  dataSet <- dataSets[ii];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                            ii, dataSet, length(dataSets)));

  chipType <- "GenomeWideSNP_6";
  combineAlleles <- FALSE;
  ds <- CnChipEffectSet$byName(dataSet, chipType=chipType, 
                           mergeStrands=TRUE, combineAlleles=combineAlleles);
  setFullNameTranslator(ds, function(names, ...) {
    gsub(crmaTags, names(crmaTags), names, fixed=TRUE);
  });

  verbose && print(verbose, ds);
  # Nothing to do?
  if (nbrOfFiles(ds) == 0) {
    verbose && cat(verbose, "No TCGA data files. Nothing to import.");
    verbose && exit(verbose);
    next;  
  }

  dsList[[ii]] <- ds;
  
  verbose && exit(verbose);
} # for (ii ...)

# Drop empty
ns <- sapply(dsList, FUN=length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for all sample names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patterns <- BiospecimenCoreResource$getBarcodePatterns();
pattern <- sprintf("^(%s).*", patterns$patient);
sampleNamesList <- lapply(dsList, FUN=function(ds) {
  names <- getNames(ds);
  sort(unique(names));
});

sampleNames <- sort(unique(unlist(sampleNamesList, use.names=FALSE)));
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1));
}
stopifnot(length(sampleNames) > 0);


# Drop files that do not contain the samples of interest
for (kk in seq(along=sampleNamesList)) {
  sampleNamesSet <- sampleNamesList[[kk]];
  ds <- dsList[[kk]];
  keep <- sapply(sampleNamesSet, FUN=function(names) {
    any(is.element(names, sampleNames));
  });
  ds <- extract(ds, keep);
  dsList[[kk]] <- ds;
} # for (kk ...)

# Drop data sets that do not contain the samples of interest
ns <- sapply(dsList, FUN=length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);

verbose && cat(verbose, "Samples of interest:");
verbose && print(verbose, sampleNames);

samplePatterns <- sprintf("^%s", sampleNames);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ii in seq(along=dsList)) {
  ds <- dsList[[ii]];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                      ii, getName(ds), length(dsList)));
  verbose && print(verbose, ds);

  # Import each data file
  verbose && enter(verbose, "Exporting TCGA data as aroma files");
  ds2 <- exportTotalAndFracB(ds, verbose=verbose);
  verbose && print(verbose, ds2);
  verbose && exit(verbose);
  
  verbose && exit(verbose);
} # for (ii ...)



############################################################################
# HISTORY:
# 2010-01-07
# o Created.
############################################################################
