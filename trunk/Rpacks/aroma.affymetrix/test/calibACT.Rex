library(aroma.affymetrix)
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- filePath("raw", expandLinks="any");
dataSets <- list.files(pattern="_(199|200)[0-9]", path=path, full.names=TRUE);
ans <- menu(choices=dataSets, title="Select data set to calibrate")
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
  path <- filePath(dataSet, chipType, expandLinks="any");
  if (!isDirectory(path))
    next;
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  calib <- AllelicCrosstalkCalibrator(ds);
  dsC <- process(calib, verbose=less(verbose));
  print(dsC);
} # for (chipType ...)

