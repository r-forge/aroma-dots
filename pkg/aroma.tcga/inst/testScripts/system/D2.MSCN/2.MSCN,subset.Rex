############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Description:
# 1. Finds all data sets with CN ratios.
# 2. Let user select by tumor type, centers, method, and data set names.
# 3. Runs MSCN on each remaining sample.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# Requires aroma.cn v0.4.5 or newer
stopifnot(!isOlderThan(aroma.cn, "0.4.5"));

figForce <- FALSE;

timeStamp <- format(Sys.time(), format="%Y-%m-%d %k:%M:%S %Z");
pid <- sprintf("%s,%s", System$getHostname(), timeStamp);

analysisTag <- "HB20091001";
analysisTag <- "HB20100105";
analysisTag <- "HB20100107";
analysisTag <- "HB20100114";

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

# MSCN requires at least to data sources
stopifnot(length(dataSets) >= 2);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  ans <- selectMenu(dataSets, selected=TRUE, title="Select two or more data sets:");
  keep <- is.element(dataSets, ans);
  dsList <- dsList[keep];
  stopifnot(length(dataSets) >= 2);
}

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
# Target loci used for fitting the normalization curve
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb,HB20090503");
verbose && print(verbose, fitUgp);

dsList <- GenericDataFileSetList(dsList);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process sample by sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=sampleNames)) {
  sampleName <- sampleNames[kk];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                   kk, sampleName, length(sampleNames)));

  dfList <- getFileList(dsList, sampleName, dropMissing=TRUE);

  # Sanity check - assert that multiple preprocessing methods not are used
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

  # Use small degrees of freedom (df=3.5; default df=5).
  # This solves the issue with Illumina MSCN data for the Ovarian
  # TCGA-09-1670 (see email with Jun Li, 2009-11-19--2010-01-14).
  df <- 3.5;
  dfTag <- sprintf("df=%.1f", df);
  mscn <- MultiSourceCopyNumberNormalization(dsListKK, fitUgp=fitUgp, 
             targetDimension=targetDimension, align="byChromosome", 
                                         tags=c("*", analysisTag, dfTag));
  mscn$.pcBandwidth <- df;
  verbose && print(verbose, mscn);
  
  dsNList <- process(mscn, verbose=verbose);
  verbose && print(verbose, dsNList);


  verbose && enter(verbose, "Doing MSCN QC plots");
  name <- "TCGA";
  tags <- getTags(mscn, collapse=",");
  path <- file.path("reports", name, tags, "QC,pairs");
  path <- Arguments$getWritablePath(path);

  allNames <- getAllNames(mscn);
  for (ii in seq(along=allNames)) {
    sampleName <- allNames[ii];
    verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", 
                                 ii, sampleName, length(allNames)));

    figName <- sprintf("%s,pairs", sampleName);
    filename <- sprintf("%s.png", figName);
    pathname <- file.path(path, filename);
    if (figForce || !isFile(pathname)) {
      devNew("png", pathname, width=640, height=640);
      par(mar=c(3.2,4.2,4.0,1.1), mgp=c(3,0.6,0), cex.axis=1.4, tcl=-0.2);
      plotSmoothedPairsOne(mscn, sampleName=sampleName, 
                                                 verbose=less(verbose, 5));
      devDone();
    }

    figName <- sprintf("%s,pairs,backtransforms", sampleName);
    filename <- sprintf("%s.png", figName);
    pathname <- file.path(path, filename);
    if (figForce || !isFile(pathname)) {
      fit <- fitOne(mscn, sampleName, verbose=less(verbose, 5));
      devNew("png", pathname, width=640, height=640);
      par(mar=c(4.6,4.2,3.0,2.1), mgp=c(2,0.6,0), cex.axis=1.4, tcl=-0.2);
      plotBacktransforms(fit);
      title(main=sampleName, outer=TRUE, line=-1, cex.sub=1.2);
      devDone();
    }

    verbose && exit(verbose);
  } # for (ii ...)
  verbose && exit(verbose);

  verbose && exit(verbose);
} # for (kk ...)



############################################################################
# HISTORY:
# 2010-01-14
# o Added automatic generation of QCs after MSCN normalization.
# o Changed the principal curve smoothness.
# 2010-01-05
# o Created from test20090929,MSCN.Rex.
############################################################################
