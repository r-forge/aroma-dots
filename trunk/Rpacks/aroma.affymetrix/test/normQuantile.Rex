library(aroma.affymetrix)
source("init.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- sapply(c("raw", "calibAllelicCT"), FUN=filePath, expandLinks="any");
dataSets <- list.files(pattern="_(199|200)[0-9]", path=path, full.names=TRUE);
ans <- menu(choices=dataSets, title="Select data set to normalize")
dataSet <- dataSets[ans];

path <- dataSet;
pattern <- "^(Mapping|BI_SNP)"
chipTypes <- list.files(pattern=pattern, path=path);

cat("----------------------------------------------------\n");
cat("Data set: ", dataSet, "\n", sep="");
cat("Path: ", path, "\n", sep="");
cat("Chip types: ", paste(chipTypes, collapse=", "), "\n", sep="");
cat("\n");

for (chipType in chipTypes) {
  if (FALSE) {
    # Reference set
    dataSetName <- "Affymetrix_2006-HapMap270,ACT";
    path <- filePath(dataSet, chipType);
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
  path <- filePath(dataSet, chipType);
  if (!isDirectory(path))
    next;

  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  norm <- QuantileNormalizer(ds, subsetToAvg=1/3, targetDistribution=td);
  dsN <- process(norm, verbose=less(verbose));
  print(dsN);
} # for (chipType ...)

