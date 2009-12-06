if (interactive()) savehistory();
library("aroma.tcga");
library("R.menu");
sourceDirectory("R/", modifiedOnly=TRUE);
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

rootPath <- "rawData/";
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
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  center <- textMenu(uCenters, value=TRUE);
} else {
  center <- "hudsonalpha.org";
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, center);
dataSets <- dataSets[centers == center];


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
  chipType <- uChipTypes[1];
  chipType <- "Human1M-Duo";
  chipType <- "HumanHap550";
}
verbose && cat(verbose, "Chip type:");
verbose && print(verbose, chipType);
keep <- sapply(chipTypeList, FUN=function(chipTypes) { 
  any(chipTypes == chipType);
});
dataSets <- dataSets[keep];


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
    "5NN" = list(
      clazz = BroadTotalCopyNumberTcgaDataFileSet,
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
} else {
  throw("Unsupported center: ", center);
}
if (interactive()) {
  typeName <- textMenu(names(typeList), value=TRUE);
} else {
  typeName <- names(typeList)[1];
  typeName <- "XandYintensity";
#  typeName <- "B_Allele_Freq";
}
verbose && cat(verbose, "Type of data: ", typeName);
type <- typeList[[typeName]];

clazz <- type$clazz;
exportFcn <- type$exportFcn;
if (is.null(exportFcn)) {
  throw("No export function specified/implemented: ", typeName);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get unit names file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Get unit names file");
if (chipType == "GenomeWideSNP_6") {
  require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
  unf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6");
} else if (chipType == "HumanHap550") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="adf,TCGA,HB20091107");
} else if (chipType == "Human1M-Duo") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="TCGA,OV,hudsonalpha.org,HB20090506");
} else if (chipType == "HG-CGH-244A") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="HB20090419");
} else if (chipType == "HG-CGH-415K") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="adf,TCGA,hms.harvard.edu,HB20090824");
} else {
  throw("Unknown TCGA chip type: ", chipType);
}
verbose && print(verbose, unf);
verbose && exit(verbose);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ii in seq(along=dataSets)) {
  dataSet <- dataSets[ii];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                            ii, dataSet, length(dataSets)));

  ds <- clazz$byName(dataSet, subdirs=chipType, paths=rootPath);
  setFullName(ds, dataSet);
  verbose && print(verbose, ds);
  # Nothing to do?
  if (nbrOfFiles(ds) == 0) {
    verbose && cat(verbose, "No TCGA data files. Nothing to import.");
    verbose && exit(verbose);
    next;  
  }

  # Import each data file
  verbose && enter(verbose, "Exporting TCGA data as aroma files");
  dsList <- exportFcn(ds, unf=unf, verbose=verbose);
  verbose && print(verbose, dsList);
  verbose && exit(verbose);
  
  verbose && exit(verbose);
} # for (ii ...)



############################################################################
# HISTORY:
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
