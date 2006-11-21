library(aroma.affymetrix);
source("init.R");

moreUnits <- 30;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("MBEI,+-", "MBEI,A+B", "RMA,+-", "RMA,A+B");
ans <- menu(choices=choices, title="Select probe-level model (PLM)")
model <- choices[ans];
modelTags <- unlist(strsplit(model, split=","));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("normQuantile");
path <- sapply(paths, FUN=filePath, expandLinks="any");
dataSets <- list.files(pattern="_(199|200)[0-9]", path=path, full.names=TRUE);
ans <- menu(choices=dataSets, title="Select data set to modelled")
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
  path <- filePath(dataSet, chipType);
  if (!isDirectory(path))
    next;

  verbose && enter(verbose, "Defining data set");
  verbose && cat(verbose, "Path: ", path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);
  verbose && exit(verbose);

  if (nbrOfArrays(ds) < 6) {
    verbose && cat(verbose, "WARNING: Too small data set");
    verbose && enter(verbose, "Appending reference data set");

    verbose && enter(verbose, "Defining reference data set");
    path <- getParent(getParent(getPath(ds)));
    refTags <- paste(getTags(ds), collapse=",");
    refFullName <- paste("Affymetrix_2006-HapMap270", refTags, sep=",");
    path <- filePath(path, refFullName, chipType);
    verbose && cat(verbose, "Path: ", path);
    dsR <- AffymetrixCelSet$fromFiles(path);
    print(dsR);
    verbose && exit(verbose);

    verbose && enter(verbose, "Appending");
    append(ds, dsR);
    verbose && cat(verbose, "Extended data set:");
    print(ds);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Defining PLM model");
  if ("MBEI" %in% modelTags) {
    plm <- MbeiCnPlm;
  } else if ("RMA" %in% modelTags) {
    plm <- RmaCnPlm;
  } else {
    throw("Unknown PLM: ", paste(modelTags, collapse=","));
  }
  mergeStrands <- !("+-" %in% modelTags);
  combineAlleles <- ("A+B" %in% modelTags);
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
    uDone <- fit(model, units=units, moreUnits=moreUnits, verbose=verbose);
    nDone <- length(uDone);
  
    verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (chr ...) 

  verbose && exit(verbose);
}

