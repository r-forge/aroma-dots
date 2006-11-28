library(aroma.affymetrix);
source("init.R");
source("ui.R");

moreUnits <- 30;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("MBEI,+-", "MBEI,A+B", "RMA,+-", "RMA,A+B");
ans <- menu(choices=choices, title="Select probe-level model (PLM)")
model <- choices[ans];
modelTags <- unlist(strsplit(model, split=","));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("normQuantile");
pattern <- "_(199|200)[0-9]";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare model class
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ("MBEI" %in% modelTags) {
  Plm <- MbeiCnPlm;
} else if ("RMA" %in% modelTags) {
  Plm <- RmaCnPlm;
} else {
  throw("Unknown PLM: ", paste(modelTags, collapse=","));
}
mergeStrands <- !("+-" %in% modelTags);
combineAlleles <- ("A+B" %in% modelTags);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");
res <- lapply(dataSets, FUN=function(ds) {
  verbose && enter(verbose, "Data set ", getName(ds));

  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  if (nbrOfArrays(ds) < 3)
    throw("Too few arrays to fit PLM too");

  verbose && enter(verbose, "Defining PLM model");
  model <- Plm(ds, mergeStrands=mergeStrands, combineAlleles=combineAlleles);
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

  model;
})
verbose && exit(verbose);


############################################################################
# HISTORY:
# 2006-11-xx
# o Recreated.
############################################################################
