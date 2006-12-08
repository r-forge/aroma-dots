library(oligo);
library(aroma.affymetrix);
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "plmData";
path <- sapply(paths, FUN=filePath, expandLinks="any");
pattern <- ",(MBEI|RMA),[+]-";
dataSets <- list.files(pattern=pattern, path=path, full.names=TRUE);
ans <- menu(choices=dataSets, title="Select data set to normalize")
dataSet <- dataSets[ans];

path <- dataSet;
pattern <- "^(Mapping|BI_SNP)"
chipTypes <- list.files(pattern=pattern, path=path);

cat("----------------------------------------------------\n");
cat("Data set: ", dataSet, "\n", sep="");
cat("Path: ", path, "\n", sep="");
cat("Chip types: ", paste(chipTypes, collapse=", "), "\n", sep="");
cat("\n");

if (length(chipTypes) > 1) {
  ans <- menu(choices=chipTypes, title="Select chip types to be considered (0 for all)");
  if (ans > 0)
    chipTypes <- chipTypes[ans];
}

for (chipType in chipTypes) {
  path <- filePath(dataSet, chipType, expandLinks="any");
  if (!isDirectory(path))
    next;

  verbose && enter(verbose, "Defining SNP chip-effect set");
  ces <- SnpChipEffectSet$fromFiles(path);
  print(ces);
  verbose && exit(verbose);

  verbose && enter(verbose, "Fitting CRLMM");
  fit <- fitCrlmm(ces, verbose=verbose);
  verbose && exit(verbose);

  next;

  verbose && enter(verbose, "Fit CRLMM");
  allChromosomes <- c(1:22, "X");
  for (chr in allChromosomes) {
    verbose && enter(verbose, "Chromosome ", chr);

    verbose && enter(verbose, "Retrieving genome information");
    cdf <- getCdf(model);
    gi <- getGenomeInformation(cdf);
    units <- getUnitIndices(gi, chromosome=chr);
    nunits <- length(units);
    verbose && exit(verbose);

    verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
    uDone <- fit(model, units=units, moreUnits=1, verbose=verbose);
    nDone <- length(uDone);
  
    verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (chr ...) 

  verbose && exit(verbose);
}

