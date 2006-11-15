library(aroma.affymetrix)
source("init.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataPath <- "raw";
#dataPath <- "calibAllelicCT";

chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
chipTypes <- c("BI_SNP");

for (chipType in chipTypes) {
  if (FALSE) {
    # Reference set
    dataSetName <- "Affymetrix_2006-HapMap270,b";
    path <- filePath(dataPath, dataSetName, chipType);
    print(path);
    ds <- AffymetrixCelSet$fromFiles(path);
    print(ds);
  
    norm <- QuantileNormalizer(ds, subsetToAvg=1/3);
    dsN <- process(norm, verbose=less(verbose));
    print(dsN);
  
    td <- getTargetDistribution(norm);
  } else {
    td <- NULL;
  }

  # Test set
  dataSetName <- "CarmichaelC_etal_2006-500k,b";
  dataSetName <- "SlaterH_etal_2004";
  dataSetName <- "Affymetrix_2006b-500k_v2";
  path <- filePath(dataPath, dataSetName, chipType);
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  norm <- QuantileNormalizer(ds, subversionTag="a", subsetToAvg=1/3, targetDistribution=td);
  dsN <- process(norm, verbose=less(verbose));
  print(dsN);
} # for (chipType ...)

