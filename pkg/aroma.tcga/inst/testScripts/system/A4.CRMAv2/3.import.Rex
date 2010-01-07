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

# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);

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

  # Replace aroma.affymetrix tags with the shorter version
  setFullNameTranslator(ds, function(names, ...) {
    gsub(crmaTags, names(crmaTags), names, fixed=TRUE);
  });

  # Replace comma-separated full names with TCGA style full names, e.g.
  # 'TCGA-06-0168,10A,01D-0237-01' to 'TCGA-06-0168-10A-01D-0237-01'.
  # Yes, this is like going back, but doing so will put all data files.
  setFullNamesTranslator(ds, function(names, ...) {
    patterns <- BiospecimenCoreResource$getBarcodePatterns();
    pattern <- sprintf("(%s),(%s),(%s-%s)", patterns$patient, 
          patterns$sampleId, patterns$portionId, patterns$plateBarcode);
    gsub(pattern, "\\1-\\5-\\8", names);
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
  names <- gsub(pattern, "\\1", names);
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

  names <- getNames(ds, translate=FALSE);
  keep <- is.element(names, sampleNames);
  ds <- extract(ds, keep);

  # BUG?  Why do we have to reset the fullnames translator here?!?
##  setFullNamesTranslator(ds, function(names, ...) {
##    patterns <- BiospecimenCoreResource$getBarcodePatterns();
##    pattern <- sprintf("(%s),(%s),(%s-%s)", patterns$patient, 
##          patterns$sampleId, patterns$portionId, patterns$plateBarcode);
##    gsub(pattern, "\\1-\\5-\\8", names);
##  });

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
