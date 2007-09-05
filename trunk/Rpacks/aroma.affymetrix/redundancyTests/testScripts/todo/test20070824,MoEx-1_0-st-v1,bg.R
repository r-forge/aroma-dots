library(aroma.affymetrix);

log <- Verbose(threshold=-4, timestamp=TRUE);

dataSetName <- "MICE_2007";
chipType <- "MoEx-1_0-st-v1";  # Not yet supported!

cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
print(cs);

bg <- GcRmaBackgroundCorrection(cs);
print(bg);
csBg <- process(bg, verbose=log);
