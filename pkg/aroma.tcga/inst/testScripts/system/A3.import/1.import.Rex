############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/import/");
#
# Command line:
#  nohupR 1.import.Rex --center=<center> --chipType=<chipType> 
#    --typeName=<typeName> [--dataSets="<dataSet1>,...,<dataSetK>"
#
# Description:
# 1. ...
# 2. ...
############################################################################
if (interactive()) savehistory();
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

  # Required arguments
  stopifnot(is.element("center", names(args)));
  stopifnot(is.element("chipType", names(args)));
  stopifnot(is.element("typeName", names(args)));

  # Optional arguments
  if (!is.null(args$dataSets)) {
    args$dataSets <- trim(strsplit(args$dataSets, split=", ")[[1]]);
  }
  verbose && cat(verbose, "Parsed command line arguments:");
  verbose && print(verbose, args);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patterns <- TcgaDccData$getDataSetPatterns();
pattern <- patterns$dataset;

rootPath <- "rawData/";
rootPath <- Arguments$getReadablePath(rootPath);
dataSets <- list.files(path=rootPath, pattern=pattern);

verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  center <- textMenu(uCenters, value=TRUE);
} else {
  center <- args$center;
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, center);
dataSets <- dataSets[centers == center];

verbose && cat(verbose, "Filtered data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipTypeList <- lapply(dataSets, FUN=function(dataSet) {
  path <- file.path(rootPath, dataSet);
  dirs <- list.files(path=path);
  dirs <- gsub("[.](lnk|LNK)", "", dirs);
  dirs;
});
uChipTypes <- sort(unique(unlist(chipTypeList, use.names=FALSE)));
if (interactive()) {
  chipType <- textMenu(uChipTypes, value=TRUE);
} else {
  chipType <- args$chipType;
}
verbose && cat(verbose, "Chip type:");
verbose && print(verbose, chipType);
keep <- sapply(chipTypeList, FUN=function(chipTypes) { 
  any(chipTypes == chipType);
});
dataSets <- dataSets[keep];

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
} else {
  if (!is.null(args$dataSets)) {
    dataSets <- args$dataSets;
  }
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select type of data to import
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (center == "broad.mit.edu") {
  typeList <- list(
    "birdseed" = list(
      clazz = BroadBirdseedGenotypeTcgaDataFileSet,
      exportFcn = exportGenotypeCallsAndConfidenceScores
    ),
    "ismpolish" = list(
      clazz = BroadAlleleSpecificThetaTcgaDataFileSet
    ),
    "CN" = list(
      clazz = BroadTotalCopyNumberTcgaDataFileSet,
      exportFcn = exportTotalAndFracB,
      extPattern = "^[^.]+[.](copynumber[.]data[.]txt)$"
    ),
    "5NN" = list(
      clazz = BroadAfter5NNTotalCopyNumberTcgaDataFileSet,
      exportFcn = exportTotalAndFracB
    ),
    "ASCN" = list(
      clazz = BroadAlleleSpecificCopyNumberTcgaDataFileSet
    )
  );
} else if (center == "hudsonalpha.org") {
  typeList <- list(
    "Genotypes" = list(
      clazz = HudsonAlphaGenotypeTcgaDataFileSet,
      exportFcn = exportGenotypeCalls
    ),
    "XandYintensity" = list(
      clazz = HudsonAlphaXYTcgaDataFileSet,
      exportFcn = exportTotalAndFracB
    ),
    "B_Allele_Freq" = list(
      clazz = HudsonAlphaBAlleleFreqTcgaDataFileSet,
      exportFcn = exportTotalAndFracB
    )
  );
} else if (center == "hms.harvard.edu") {
  typeList <- list(
    "TCN" = list(
      clazz = HarvardTotalCopyNumberTcgaDataFileSet,
      exportFcn = exportTotal
    )
  );
} else if (center == "mskcc.org") {
  typeList <- list(
    "TCN" = list(
      clazz = MskccTotalCopyNumberTcgaDataFileSet,
      exportFcn = exportTotal
    )
  );
} else {
  throw("Unsupported center: ", center);
}
if (interactive()) {
  typeName <- textMenu(names(typeList), value=TRUE);
} else {
  typeName <- args$typeName;
}
verbose && cat(verbose, "Type of data: ", typeName);
type <- typeList[[typeName]];

clazz <- type$clazz;
exportFcn <- type$exportFcn;
if (is.null(exportFcn)) {
  throw("No export function specified/implemented: ", typeName);
}
extPattern <- type$extPattern;


verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup list of data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list();
for (ii in seq(along=dataSets)) {
  dataSet <- dataSets[ii];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                            ii, dataSet, length(dataSets)));

  if (is.null(extPattern)) {
    ds <- clazz$byName(dataSet, subdirs=chipType, paths=rootPath);
  } else {
    ds <- clazz$byName(dataSet, subdirs=chipType, paths=rootPath, 
                                                  pattern=extPattern);
  }
  setFullName(ds, dataSet);
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
# Scan for all data files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filenamesList <- lapply(dsList, FUN=function(ds) {
  sapply(ds, getFilename);
});

filenames <- sort(unique(unlist(filenamesList, use.names=FALSE)));
if (interactive()) {
  filenames <- selectMenu(filenames, selected=TRUE);
}
stopifnot(length(filenames) > 0);


# Drop files that do not contain the samples of interest
for (kk in seq(along=filenamesList)) {
  filenamesSet <- filenamesList[[kk]];
  ds <- dsList[[kk]];
  keep <- sapply(filenamesSet, FUN=function(names) {
    any(is.element(names, filenames));
  });
  ds <- extract(ds, keep);
  dsList[[kk]] <- ds;
} # for (kk ...)

# Drop data sets that do not contain the samples of interest
ns <- sapply(dsList, FUN=length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for all sample names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patterns <- BiospecimenCoreResource$getBarcodePatterns();
pattern <- sprintf("^(%s).*", patterns$patient);
sampleNamesList <- lapply(dsList, FUN=function(ds) {
  lapply(ds, function(df) {
    names <- getColumnNames(df);
    names <- grep(pattern, names, value=TRUE);
    names <- gsub(pattern, "\\1", names);
    sort(unique(names));
  });
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
# Get unit names file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Get unit names file");
if (chipType == "GenomeWideSNP_6") {
  require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
  unf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6");
} else if (chipType == "HumanHap550") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="TCGA,HB20100107");
} else if (chipType == "Human1M-Duo") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="TCGA,OV,hudsonalpha.org,HB20090506");
} else if (chipType == "HG-CGH-244A") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="HB20090419");
} else if (chipType == "HG-CGH-415K") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="adf,TCGA,hms.harvard.edu,HB20090824");
} else if (chipType == "HG-CGH-1x1M") {
 unf <- TextUnitNamesFile$byChipType(chipType, tags="HB20090419");
} else {
  throw("Unknown TCGA chip type: ", chipType);
}
verbose && print(verbose, unf);
verbose && exit(verbose);


# Maximum number of unknown unit names that are accepted (will be skipped)
maxNbrOfUnknownUnitNames <- 2;

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
  dsListOut <- exportFcn(ds, unf=unf, samplePatterns=samplePatterns,
       maxNbrOfUnknownUnitNames=maxNbrOfUnknownUnitNames, verbose=verbose);
  verbose && print(verbose, dsListOut);
  verbose && exit(verbose);
  
  verbose && exit(verbose);
} # for (ii ...)



############################################################################
# HISTORY:
# 2010-01-06
# o Added the option to select a subset of sample names too.
# 2010-01-03
# o Added support to import total CNs from hms.harvard.edu.
# 2009-12-06
# o Added support to import Genotypes.txt from hudsonalpha.org.
# 2009-11-06
# o Updated to also import XandYintensity.txt from hudsonalpha.org.
# 2009-11-02
# o Updated to also import B_Allele_Freq.txt from hudsonalpha.org.
# o Added support for importing broad.mit.edu Birdseed genotype calls.
# o Added support for importing hudsonalpha.org BAFs.
# 2009-09-24
# o Added verbose output.
# 2009-xx-xx
# o Created.
############################################################################
