library(aroma.affymetrix);
#source("../aroma.affymetrix/R/AllelicCrosstalkCalibration.R");
library(R.graphics);
log <- Verbose(threshold=-5, timestamp=TRUE);


imgFormat <- "screen";
imgFormat <- "eps";

fig <- 1;

name <- "Affymetrix_2006-HapMap270.CEU.founders";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");
chipType <- chipTypes[1];

if (!exists("cs")) {
  cs <- AffymetrixCelSet$fromName(name, chipType=chipType);
  cs <- extract(cs, 1:2);
  print(cs);
}


# Identify units and cell on chr X and their complements
if (!exists("acc")) {
  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  xUnits <- getUnitsOnChromosome(gi, 23);
  nonXUnits <- setdiff(seq_len(nbrOfUnits(cdf)), xUnits);
  xCells <- unlist(getCellIndices(cdf, units=xUnits), use.names=FALSE);
  nonXCells <- setdiff(seq_len(nbrOfCells(cdf)), xCells);
  
  acc <- AllelicCrosstalkCalibration(cs, subsetToAvg=nonXCells, tags=c("*","-X"));
#  acc <- AllelicCrosstalkCalibration(cs);
  rm(csC);
  print(acc);
}

if (!exists("csC")) {
  csC <- process(acc, verbose=verbose);
  print(csC);
}

accLinesAllelicCrosstalk <- function(a, B, max=1e5, lwd=2, col="black", lty="31", ...) {
str(a);
str(B);
  bA <- B[1,2]/B[1,1]; 
  bB <- B[2,1]/B[2,2]; 
  lines(x=a[1]+c(0,1)*max, y=a[2]+c(0,1)*max*bA, col="white", lwd=1.7*lwd, ...);
  lines(x=a[1]+c(0,1)*max*bB, y=a[2]+c(0,1)*max, col="white", lwd=1.7*lwd, ...);
  lines(x=a[1]+c(0,1)*max, y=a[2]+c(0,1)*max*bA, lwd=lwd, col=col, lty=lty, ...);
  lines(x=a[1]+c(0,1)*max*bB, y=a[2]+c(0,1)*max, lwd=lwd, col=col, lty=lty, ...);
}

fig <- 1;

sampleName <- "NA06985";
#sampleName <- getNames(cs)[2];
basepair <- "AT";
xlim <- c(-800,8500);
scales <- c(before=1, after=0.15);
colramp <- function(n) { grey(seq(from=1, to=0, length=n)); }
transformation <- function(x) x^0.33;
#transformation <- function(x) x^0.1;

for (what in c("before", "after")) {
  if (!Device$isOpen(fig <- fig + 1)) {
    Device$set(fig);
  
    useColorsStr <- "";
    imgName <- sprintf("%s,%s,%s", sampleName, basepair, what);
    filename <- sprintf("%s%s.eps", imgName, useColorsStr);
    eps(filename, width=6, height=6);
    par(mar=c(4.6,5.4,1.8,1.8)+0.1);
    par(cex=0.7, mgp=c(3.6,1.1,0));
    par(cex.main=2, cex.axis=2, cex.lab=2, font.lab=2);
#    par(cex.main=2, cex.axis=3, cex.lab=3, font.lab=2);
  
    plotBasepair(acc, array=1, basepair="AT", xlim=xlim, 
                 transformation=transformation,
                 colramp=colramp, scale=scales[what],
                 linesFcn=accLinesAllelicCrosstalk, lcol="black",
                 what=what, bandwidth=50,
                 verbose=log);
  
    dev.off();
  }
} # for (what ...)
