library(aroma.affymetrix)
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetName <- "CarmichaelC_etal_2006-500k";
dataSetName <- "Affymetrix_2006-HapMap270";
dataSetName <- "SlaterH_etal_2006";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

#dataSetName <- "SlaterH_etal_2004";
#chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

for (chipType in chipTypes) {
  path <- filePath("raw", dataSetName, chipType);
  if (!isDirectory(path))
    next;
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  calib <- AllelicCrosstalkCalibrator(ds);
  dsC <- process(calib, verbose=less(verbose));
  print(dsC);
} # for (chipType ...)

