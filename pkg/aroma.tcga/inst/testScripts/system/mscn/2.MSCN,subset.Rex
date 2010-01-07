############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/normalize/");
#
# Description:
# Applies MSCN normalization to a set of samples, which are specified by
# loaded data fram of pre-scanned samples.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

timeStamp <- format(Sys.time(), format="%Y-%m-%d %k:%M:%S %Z");
pid <- sprintf("%s,%s", System$getHostname(), timeStamp);

analysisTag <- "HB20091001";
analysisTag <- "HB20100105";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all tumor-normal CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "rawCnData";
rootPath <- Arguments$getReadablePath(rootPath, mustExist=TRUE);
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,.*", TcgaDccData$getDataSetPatterns()$dataset);
types <- gsub(pattern, "\\2", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select type:");
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select preprocessing types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,", TcgaDccData$getDataSetPatterns()$dataset);
types <- gsub(pattern, "", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select type:");
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select samples to process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get all data sets
dsList <- lapply(dataSets, FUN=function(dataSet) {
  AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
});

dsList <- GenericDataFileSetList(dsList, allowDuplicates=FALSE);
sampleNames <- getNames(dsList);

# User selection
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1));
}

dsList <- extract(dsList, sampleNames);
sampleNames <- getNames(dsList);

dataSets <- unname(sapply(getSets(dsList), FUN=getFullName));
verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Target loci used for fitting the normalization curve
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb,HB20090503");
verbose && print(verbose, fitUgp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process sample by sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=sampleNames)) {
  sampleName <- sampleNames[kk];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                   kk, sampleName, length(sampleNames)));

  dfList <- getFileList(dsList, sampleName, dropMissing=TRUE);

  # Sanity check - assert that multiple preprocessing methods are used
  pattern <- sprintf("(%s).*", TcgaDccData$getDataSetPatterns()$dataset);
  dataSetsKK <- gsub(pattern, "\\1", names(dfList));
  dups <- which(duplicated(dataSetsKK));
  ndups <- length(dups);
  stopifnot(ndups == 0);

  nbrOfSources <- length(dfList);
  verbose && cat(verbose, "Number of sites that processed this sample: ", nbrOfSources);
  if (nbrOfSources == 1) {
    verbose && cat(verbose, "Nothing to do. This sample was only process by one site: ", sampleName);
    verbose && exit(verbose);
    next;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Selecting target dimension
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Identifying prioritized target dimension");
  # Identify (center, chipType):s
  pattern <- sprintf("^%s,.*", TcgaDccData$getDataSetPatterns()$dataset);
  centers <- gsub(pattern, "\\1", names(dfList));
  chipTypes <- unname(sapply(dfList, FUN=getChipType));
  data <- data.frame(center=centers, chipType=chipTypes);
  verbose && print(verbose, data);

  # Identify target dimension
  keyOrder <- c("mskcc.org/HG-CGH-1x1M", "broad.mit.edu/GenomeWideSNP_6", 
                "hms.harvard.edu/HG-CGH-244A", "hudsonalpha.org/Human1M-Duo",
                "hms.harvard.edu/HG-CGH-415K", "mskcc.org/HG-CGH-244A");
  keys <- paste(data$center, data$chipType, sep="/");
  idxs <- match(keyOrder, keys);
  idxs <- na.omit(idxs);
  targetDimension <- idxs[1];
  stopifnot(is.finite(targetDimension));
  verbose && cat(verbose, "Target (center, chipType):");
  verbose && print(verbose, data[targetDimension,,drop=FALSE]);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Multi-source CN normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dsListKK <- lapply(dfList, FUN=function(df) {
    AromaUnitTotalCnBinarySet(list(df));
  });
  mscn <- MultiSourceCopyNumberNormalization(dsListKK, fitUgp=fitUgp, 
             targetDimension=targetDimension, align="byChromosome", 
                                         tags=c("*", analysisTag));
  verbose && print(verbose, mscn);
  
  dsNList <- process(mscn, verbose=verbose);
  verbose && print(verbose, dsNList);

  verbose && exit(verbose);
} # for (kk ...)



############################################################################
# HISTORY:
# 2010-01-05
# o Created from test20090929,MSCN.Rex.
############################################################################
