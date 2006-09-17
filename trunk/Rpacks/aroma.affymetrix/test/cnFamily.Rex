source("init.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- findCelSet(chipType="Xba", minCount=5);

if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(path);
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalize
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("qn", mode="list")) {
  verbose && enter(verbose, "Normalizing all data sets");
  tmp <- raw;
  for (name in names(raw)) {
    verbose && enter(verbose, "Data set: ", name);
    ds <- raw[[name]];
    dsN <- normalizeQuantile(ds, verbose=verbose);
    print(dsN);
    tmp[[name]] <- dsN;
    verbose && exit(verbose);
  }
  qn <- tmp;
  verbose && exit(verbose);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Estimate copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chr <- "8";
for (ds in qn) {
  verbose && enter(verbose, "Data set: ", getName(ds));
  model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);  
  verbose && print(verbose, model);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);
  units <- getUnitIndices(gi, chromosome=chr);
  nunits <- length(units);
  verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
  uDone <- fit(model, units=units, verbose=verbose);
  nDone <- length(uDone);
  verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
  verbose && exit(verbose);
  verbose && exit(verbose);
}
