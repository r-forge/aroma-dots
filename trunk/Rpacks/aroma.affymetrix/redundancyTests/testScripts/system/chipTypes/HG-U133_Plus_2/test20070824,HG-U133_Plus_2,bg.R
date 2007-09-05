library(aroma.affymetrix);

log <- Verbose(threshold=-4, timestamp=TRUE);

dataSetName <- "heart_brain";
chipType <- "HG-U133_Plus_2";

cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
print(cs);

bg <- GcRmaBackgroundCorrection(cs);
print(bg);
csBg <- process(bg, verbose=log);
