############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/normalize/");
#
# Description:
# Normalized (tumor,normal) CN ratios using principal-curve normalization.
#  1. ...
#  2. ...
############################################################################
if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getTypeTag <- function(this, ...) {
  gsub("[A-Z]$", "", getTags(this)[1]);
}

tcgaFNT01 <- function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User selects: data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);

verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

if (interactive()) {
  dataSet <- textMenu(dataSets, value=TRUE); 
} else {
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.9.5.0,ov,ismpolish";
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.11.5.0,ov,ismpolish";
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0,ov,ismpolish";
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.13.6.0,ov,ismpolish";
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.14.7.0,ov,ismpolish";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0,ov";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.2.1.0,ov";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.3.0.0,ov";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.4.0.0,ov";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.5.0.0,ov";
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.6.0.0,ov";
}

verbose && cat(verbose, "Selected data set:");
verbose && print(verbose, dataSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup total locus data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- NULL;
dsR <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*", paths=rootPath);
setFullNamesTranslator(dsR, tcgaFNT01);

verbose && cat(verbose, "Data set:");
verbose && print(verbose, dsR);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract all tumors
dsT <- extractByTcgaType(dsR, pattern="^01");

# Extract all normals
dsN <- extractByTcgaType(dsR, pattern="^1[01]");

# Keep only samples for which there is both a tumor and a normal
sampleNames <- intersect(getNames(dsT), getNames(dsN));
verbose && cat(verbose, "Identified tumor-normal samples:");
verbose && print(verbose, sampleNames);

dsT <- extract(dsT, sampleNames);
dsN <- extract(dsN, sampleNames);

# Drop duplicated samples
dsT <- extract(dsT, !duplicated(getNames(dsT)));
dsN <- extract(dsN, !duplicated(getNames(dsN)));

verbose && cat(verbose, "Data set of tumors:");
verbose && print(verbose, dsT);
verbose && cat(verbose, "Data set of normals:");
verbose && print(verbose, dsN);

# Sanity check
stopifnot(identical(getNames(dsT), getNames(dsN)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Principal curve normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  choises <- c(TRUE, FALSE);
  onLog2Scale <- textMenu(choises, title="Fit on log2 scale?", value=TRUE); 
} else {
  onLog2Scale <- TRUE;
}

setFullNamesTranslator(dsT, NULL);
setFullNamesTranslator(dsN, NULL);
pcn <- PrincipalCurveNormalization(dsT, dsN);
# Fit on the log2 scale?
if (onLog2Scale) {
  pcn$.transformFcn <- function(x) log2(x);
  pcn$.untransformFcn <- function(y) 2^y;
  setTags(pcn, c("*", "log2"));
}
verbose && print(verbose, pcn);

dsN <- process(pcn, verbose=verbose);

verbose && cat(verbose, "Normalized data set:");
verbose && print(verbose, dsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chr <- 7;
ugp <- getAromaUgpFile(dsR);
units <- getUnitsOnChromosome(ugp, chromosome=chr);
pos <- getPositions(ugp, units=units);

lim <- c(0,6);
Alim <- c(-5,3);
Mlim <- c(-3,3);

fullnames <- getFullNames(dsR);
fullnames <- fullnames[idxs[,1]];
dsList <- list(dsR, dsN);
subplots(3*length(dsList), ncol=3, byrow=TRUE);
for (ds in dsList) {
  dsPair <- extract(ds, indexOf(ds, fullnames));
  theta <- extractMatrix(dsPair, units=units);

  plot(theta[,2:1], pch=".", xlim=lim, ylim=lim);
  abline(a=0, b=1, lwd=2, col="#999999");
  subset <- sample(1:nrow(theta), size=nrow(theta)/8);
  fit <- fitPrincipalCurve(theta[subset,]);
  xy <- fit$s[,2:1];
  xy <- xy[order(xy[,1]),];
  lines(xy, col="red", lwd=2);

  plotMvsA(theta, pch=".", Alim=Alim);
  abline(h=0, lwd=2, col="#999999");

  cn <- RawCopyNumbers(cn=log2(theta[,2]/theta[,1]), x=pos);
  plot(cn, pch=".");
  stext(side=3, pos=1, sprintf("Chr %d", chr));
  cnS <- binnedSmoothing(cn, by=500e3);
  points(cnS, pch=".", cex=3, col="red");
} # for (ds ...)


############################################################################
# HISTORY:
# 2010-01-05
# o Updated script so that one can fit on the log2-scale (default).
# o Updated code for identifying tumor-normal pairs.
# o Added instructions to the header.
# 2009-07-15
# o Created.
############################################################################
