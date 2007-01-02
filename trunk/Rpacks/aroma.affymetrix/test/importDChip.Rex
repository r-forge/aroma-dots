library(aroma.affymetrix);
source("init.R");
source("ui.R");

rootPath <- "dChipData";

mergeStrands <- TRUE;
combineAlleles <- TRUE;

cdf <- AffymetrixCdfFile$fromChipType("Mapping250K_Nsp");

path <- filePath(rootPath, "SinclairA_etal_2006,QN,MBEI,PM-MM,A+B", getChipType(cdf));

filenames <- list.files(path=path, pattern=".* expression[.]xls$");
filename <- filenames[1];

ces <- CnChipEffectSet$importFromDChip(filename, path=path, combineAlleles=TRUE, cdf=cdf, skip=FALSE, verbose=verbose);


############################################################################
# HISTORY:
# 2007-01-02
# o Created.
############################################################################
