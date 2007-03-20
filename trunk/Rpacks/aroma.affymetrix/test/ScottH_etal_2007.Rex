library(aroma.affymetrix);
#source("init.R");

dataSetName <- "AGRF_2007a";
tags <- "QN,RMA,A+B,FLN";
#tags <- "QN,RMA,+-,A+B,F+R,g,FLN";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

cesList <- list();
for (chipType in chipTypes) {
  ces <- CnChipEffectSet$fromName(dataSetName, tags=tags, chipType=chipType);
  print(ces);
  cesList[[chipType]] <- ces;
} # for (chipType ...)


# Create a log file
glad <- GladModel(cesList);
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

