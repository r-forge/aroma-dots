# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
savehistory();
library(aroma.affymetrix);
library(R.graphics);
log <- Verbose(threshold=-5, timestamp=TRUE);

name <- "Affymetrix_2006-HapMap270.CEU.founders";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");
chipType <- chipTypes[1];

if (!exists("cs")) {
  cs <- AffymetrixCelSet$fromName(name, chipType=chipType);
  n23 <- as.integer(sapply(cs, FUN=getAttribute, "n23"));
}
print(cs);

# Identify units and cell on chr X and their complements
cdf <- getCdf(cs);
gi <- getGenomeInformation(cdf);
xUnits <- getUnitsOnChromosome(gi, 23);
nonXUnits <- setdiff(seq_len(nbrOfUnits(cdf)), xUnits);
xCells <- unlist(getCellIndices(cdf, units=xUnits), use.names=FALSE);
nonXCells <- setdiff(seq_len(nbrOfCells(cdf)), xCells);

if (TRUE) {
  acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=nonXCells, tags=c("*","-X"));
} else {
  acc <- AllelicCrosstalkCalibration(cs);
}
print(acc);
csC <- process(acc, verbose=log);
print(csC);

sampleNames <- getNames(cs);
poorSamples <- match(c("NA12145"), sampleNames);

col <- "black";
xlim <- c(0, 16);
ylim <- c(0, 0.4);

fig <- 1;
for (cs0 in list(cs, csC)) {
  if (!Device$isOpen(fig <- fig + 1)) {
    Device$set(fig);

    imgName <- sprintf("%s,plotDensity", getFullName(cs0));
    imgName <- gsub("[.]", ",", imgName);
    filename <- sprintf("%s.eps", imgName);
    eps(filename, width=5, height=3);
  
    par(cex.axis=1.2, cex.lab=1.2, font.lab=2);
    par(cex=0.7, mgp=c(2.1,0.8,0));
    par(mar=c(3.5,3.2,1,1)+0.1);
    
    plotDensity(cs0, types="pm", subset=1/8, n=512/4, col=col, lty=1, xlim=xlim, ylim=ylim, annotate=FALSE, verbose=log)
    
    dev.off();
  }
}
