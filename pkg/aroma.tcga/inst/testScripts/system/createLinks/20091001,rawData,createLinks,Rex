############################################################################
# Author: Henrik Bengtsson, hb@stat.berkeley.edu
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/createLinks/");
#
# Description:
# The purpose of this script is to create data sets under rawData/ that
# links (via file system links) to the corresponding tcgaData/ directory.
# For example rawData/mskcc.org_OV.CGH-1x1M_G4447A.18.1.0/HG-CGH-1x1M ->
# tcgaData/anonymous/tumor/ov/cgcc/mskcc.org/cgh-1x1m_g4447a/cna/
# mskcc.org_OV.CGH-1x1M_G4447A.18.1.0/.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for all TCGA data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find all TCGA data sets
srcPaths <- TcgaDccData$findDataSets();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extra chip types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSets <- basename(srcPaths);
pattern <- "^([^_]+)_([^.]+)[.]([^.]+)[.]([0-9]+[.][0-9]+[.][0-9])$";
chipTypes <- gsub(pattern, "\\3", dataSets);

knownChipTypes <- c(
  "HG-CGH-244A"        = "HG-CGH-244A",
  "CGH-1x1M_G4447A"    = "HG-CGH-1x1M",
  "HG-CGH-415K_G4124A" = "HG-CGH-415K",
  "Genome_Wide_SNP_6"  = "GenomeWideSNP_6",
  "HumanHap550"        = "HumanHap550",
  "Human1MDuo"         = "Human1M-Duo"
);

# Translate
chipTypes <- knownChipTypes[chipTypes];

# Sanity check
stopifnot(!any(is.na(chipTypes)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create data set directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "rawData";
rootPath <- Arguments$getWritablePath(rootPath);
destPaths <- file.path(rootPath, dataSets);
sapply(destPaths, FUN=mkdirs);
# Sanity check
sapply(destPaths, FUN=function(path) Arguments$getWritablePath(path));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create (missing) chip type directories linking to TCGA directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
linkPaths <- file.path(destPaths, chipTypes);
for (kk in seq(along=linkPaths)) {
  link <- linkPaths[kk];
  target <- srcPaths[kk];
  pathname <- Arguments$getReadablePathname(link, mustExist=FALSE);
  if (!file.exists(pathname)) {
    createLink(link, target=target);
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate by setting up data file sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=dataSets)) {
  dataSet <- dataSets[kk];
  chipType <- chipTypes[kk];
  ds <- GenericDataFileSet$byName(dataSet, subdirs=chipType, paths=rootPath);
  print(ds);
}


############################################################################
# HISTORY:
# 2010-01-02
# o Added a description.
# 2009-10-26
# o Removed requirement on affxparser by using TcgaDccData$findDataSets().
# 2009-10-01
# o Created.
############################################################################
