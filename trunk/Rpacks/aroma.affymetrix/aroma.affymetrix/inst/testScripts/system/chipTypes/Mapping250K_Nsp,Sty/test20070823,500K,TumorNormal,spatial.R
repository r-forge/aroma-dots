library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "Affymetrix_2006-TumorNormal";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");
#chipTypes <- chipTypes[2];

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRawList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$fromName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
  csRawList[[chipType]] <- cs;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(csRawList);
setColorMaps(ae, "log2,yellow");
if (interactive())
  process(ae, verbose=log);
