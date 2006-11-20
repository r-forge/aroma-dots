# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("MBEI,+-", "MBEI,A+B", "RMA,+-", "RMA,A+B");
ans <- menu(choices=choices, title="Select probe-level model (PLM)")
model <- choices[ans];
tags <- unlist(strsplit(model, split=","));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("normQuantile");
path <- sapply(paths, FUN=filePath, expandLinks="any");
dataSets <- list.files(pattern="_(199|200)[0-9]", path=path, full.names=TRUE);
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


for (chipType in chipTypes) {
  path <- filePath(dataSet, chipType);
  if (!isDirectory(path))
    next;

  verbose && enter(verbose, "Defining data set");
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);
  verbose && exit(verbose);

  verbose && enter(verbose, "Defining PLM model");
  if ("MBEI" %in% tags) {
    plm <- MbeiCnPlm;
  } else if ("RMA" %in% tags) {
    plm <- RmaCnPlm;
  } else {
    throw("Unknown PLM: ", paste(tags, collapse=","));
  }
  mergeStrands <- !("+-" %in% tags);
  combineAlleles <- ("A+B" %in% tags);
  model <- plm(ds, mergeStrands=mergeStrands, combineAlleles=combineAlleles);
  print(model);
  verbose && exit(verbose);

  verbose && enter(verbose, "Fit model");
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

