library(aroma.affymetrix);
source("init.R");
source("ui.R");

rootPath <- "dChipData";

mergeStrands <- TRUE;
combineAlleles <- TRUE;

dataSet <- list();
for (chipType in c("Mapping250K_Nsp", "Mapping250K_Sty")) {
  cdf <- AffymetrixCdfFile$fromChipType(chipType);

  path <- filePath(rootPath, "SinclairA_etal_2006,QN,MBEI,PM-MM,A+B", getChipType(cdf));

  filenames <- list.files(path=path, pattern=".* expression[.]xls$");
  filename <- filenames[1];

  ces <- CnChipEffectSet$importFromDChip(filename, path=path, combineAlleles=TRUE, cdf=cdf, skip=FALSE, verbose=verbose);

  dataSet[[chipType]] <- ces;
}


############################################################################
# HISTORY:
# 2007-01-02
# o Created.
############################################################################
