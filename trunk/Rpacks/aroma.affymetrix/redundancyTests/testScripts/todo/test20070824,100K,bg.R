library(aroma.affymetrix);
source("../aroma.affymetrix/R/AffymetrixCdfFile.computeAffinities.R");

log <- Verbose(threshold=-4, timestamp=TRUE);

dataSetName <- "Affymetrix_2004-100k_trios,testSet";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");


csRawList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
  keep <- 1:3;
  cs <- extract(cs, keep);
  sampleNames <- sampleNames[keep];
  print(cs);
  stopifnot(identical(getNames(cs), sampleNames));
  csRawList[[chipType]] <- cs;
}


csRmaBgList <- list();
for (chipType in chipTypes) {
  cs <- csRawList[[chipType]];
  bg <- RmaBackgroundCorrection(cs);
  print(bg);
  csRmaBgList[[chipType]] <- process(bg, verbose=log);
}


csGcRmaBgList <- list();
for (chipType in chipTypes) {
  cs <- csRawList[[chipType]];
  bg <- GcRmaBackgroundCorrection(cs);
  print(bg);
  csGcRmaBgList[[chipType]] <- process(bg, verbose=log);
}


