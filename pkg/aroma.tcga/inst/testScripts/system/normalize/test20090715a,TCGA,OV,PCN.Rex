if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tcgaFNT01 <- function(names, ...) {
  pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
  gsub(pattern, "\\1,\\2,\\3", names);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# User selects: data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);

dataSets <- sort(list.files(pattern=",ov", path=rootPath));

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


# Extract tuples of arrays for each unique sample
sampleNames <- unique(getNames(dsR));
idxs <- indexOf(dsR, sampleNames);
if (!is.matrix(idxs)) {
  # Drop non-pairs
  ns <- sapply(idxs, length);
  print(table(ns));
  isPair <- (ns == 2);
  idxs <- indexOf(dsR, sampleNames[isPair]);
  stopifnot(is.matrix(idxs));
}
stopifnot(nrow(idxs) == 2);
dimnames(idxs) <- names(idxs) <- NULL;
colnames(idxs) <- getNames(dsR)[idxs[1,]];

dsR <- extract(dsR, idxs);
idxs <- matrix(seq(dsR), nrow=2);

# Identify sample type, e.g. tumor or normal etc.
types <- sapply(dsR, function(df) getTags(df)[1]);
types <- gsub("[A-Z]$", "", types);
dim(types) <- dim(idxs);
uTypes <- unique(as.vector(types));

verbose && cat(verbose, "Unique data type (01,10,11,...):");
verbose && print(verbose, uTypes);

# Sort each column according to sample type
o <- apply(types, MARGIN=2, FUN=order);
for (kk in seq(length=ncol(o))) {
  types[,kk] <- types[o[,kk],kk];
  idxs[,kk] <- idxs[o[,kk],kk];
}
stopifnot(all(is.finite(idxs)));

dsT <- extract(dsR, idxs[1,]);
dsN <- extract(dsR, idxs[2,]);

verbose && cat(verbose, "Tumor data set:");
verbose && print(verbose, dsT);
verbose && cat(verbose, "Normal data set:");
verbose && print(verbose, dsN);

# Sanity check
stopifnot(identical(getNames(dsT), getNames(dsN)));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setFullNamesTranslator(dsT, NULL);
setFullNamesTranslator(dsN, NULL);
pcn <- PrincipalCurveNormalization(dsT, dsN);
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
# 2009-07-15
# o Created.
############################################################################
