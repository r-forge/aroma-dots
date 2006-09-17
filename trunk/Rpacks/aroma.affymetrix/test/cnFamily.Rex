source("init.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  path <- findCelSet(chipType="Xba", minCount=5);
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
