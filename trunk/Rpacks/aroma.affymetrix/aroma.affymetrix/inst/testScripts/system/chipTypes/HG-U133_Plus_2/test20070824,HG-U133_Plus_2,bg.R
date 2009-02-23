library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSetName <- "Affymetrix-HeartBrain";
chipType <- "HG-U133_Plus_2";

cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
print(cs);

bg <- GcRmaBackgroundCorrection(cs);
print(bg);
csBg <- process(bg, verbose=log);
