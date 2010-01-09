############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Command line:
#  nohupR 1.PairedTumorNormalBAFs.Rex
#    --types=<types> [--dataSets="<dataSet1>,...,<dataSetK>"]
#
# Description:
# 1. Identifies all tumor-normal pairs.
# 2. For each pair, generates a PNG with tumor vs normal BAFs
#    for autosomal chromosomes.
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



reportPath <- Arguments$getWritablePath("reports");

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
# Keep only SNP arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tcgaChipTypes <- c("Genome_Wide_SNP_6", "HumanHap550", "Human1MDuo");
pattern <- sprintf("(%s)", paste(tcgaChipTypes, collapse="|"));
dataSets <- grep(pattern, dataSets, value=TRUE);
verbose && cat(verbose, "Available SNP data sets:");
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
  dsR <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", 
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
} else {
  if (!is.null(args$sampleNames)) {
    sampleNames <- intersect(sampleNames, args$sampleNames);
  }
}

verbose && print(verbose, dsPairList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generating QC plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qcTag <- "QC,PairedTumorNormalBAFs";
for (kk in seq(along=dsPairList)) {
  dsPair <- dsPairList[[kk]];
  dsT <- dsPair$tumor;
  dsN <- dsPair$normal;
  stopifnot(getNames(dsT) == getNames(dsN));

  dataSetName <- getFullName(dsT);

  ugp <- getAromaUgpFile(dsT);
  chipType <- getChipType(ugp);
  figPath <- file.path(reportPath, dataSetName, qcTag, chipType);
  figPath <- Arguments$getWritablePath(figPath);

  # Identify units to plot
  units <- getUnitsOnChromosomes(ugp, chromosomes=1:22);
  
  for (ii in seq(dsT)) {
    dfT <- getFile(dsT, ii);
    sampleName <- getName(dfT);
    figName <- sprintf("%s,BAFvsBAF", sampleName);
    filename <- sprintf("%s.png", figName);
    pathname <- file.path(figPath, filename);
    if (!isFile(pathname)) {
      dfN <- getFile(dsN, ii);
  
      # Extract data
      yT <- dfT[units,1,drop=TRUE];
      yN <- dfN[units,1,drop=TRUE];
      Y <- matrix(c(yN, yT), ncol=2, byrow=FALSE);
      names <- c(getFullName(dfN),getFullName(dfT));
      names <- gsub(",fracB", "", names, fixed=TRUE);
      names <- gsub(",", "-", names, fixed=TRUE);
      colnames(Y) <- names;

      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
  
      # Plot
      devNew("png", pathnameT, width=640, height=640);
      par(mar=c(3.5,3.2,1.2,1.2)+0.1, mgp=c(1.9,0.6,0), 
          cex.axis=1.5, cex.lab=1.5);
      plot(Y, pch=".", xlim=c(0,1), ylim=c(0,1), axes=FALSE);
      axis(side=1, at=c(0,1/2,1));
      axis(side=2, at=c(0,1/2,1));
      box();
      stext(side=3, pos=0, dataSetName);
      stext(side=4, pos=0, sprintf("Chr01-22 (n=%d)", nrow(Y)));
      stext(side=4, pos=1, chipType);
      devDone();

      file.rename(pathnameT, pathname);
      if (isFile(pathnameT) || !isFile(pathname)) {
        throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
      }
    }
  } # for (ii ...)
} # for (kk ...)


############################################################################
# HISTORY:
# 2010-01-08
# o Created.
############################################################################
