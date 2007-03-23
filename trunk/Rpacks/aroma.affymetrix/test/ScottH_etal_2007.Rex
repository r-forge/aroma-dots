library(aroma.affymetrix);
#source("init.R");

verbose <- Arguments$getVerbose(-5);

dataSetName <- "AGRF_2007a";
tags <- "QN,RMA,A+B,FLN";
#tags <- "QN,RMA,+-,A+B,F+R,g,FLN";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

# Sample annotation data (attributes)
saf <- SampleAnnotationSet$fromPath("rawData");

cesList <- refList <- list();
for (chipType in chipTypes) {
  ces <- CnChipEffectSet$fromName(dataSetName, tags=tags, chipType=chipType);
  # Added sample annotation data as attributes
  setAttributesBy(ces, saf);
  print(ces);

  # Get baseline
  ref <- calculateBaseline(ces, defaultPloidy=2, verbose=verbose);
  print(ref);

  cesList[[chipType]] <- ces;
  refList[[chipType]] <- ref;
} # for (chipType ...)


# Create a log file
glad <- GladModel(cesList, refList);
print(glad);

ce <- ChromosomeExplorer(glad);
setAlias(ce, "ScottH_etal_2007");
setArrays(ce, c("1-288", "2-437", "3-575", "414-VA", "606-VA", "617-VA", "609-PW", "614-NB"));
print(ce);

logFile <- sprintf("%s.log", getName(ce));
log <- Verbose(logFile, threshold=-4, removeFile=FALSE);
timestampOn(log);

print(log, glad);
print(log, ce);
process(ce, verbose=log);

