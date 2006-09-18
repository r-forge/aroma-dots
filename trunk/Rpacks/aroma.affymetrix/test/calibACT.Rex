# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("dataPath"))
  dataPath <- findCelSet(chipType="Xba", minCount=5);
#dataPath <- findCelSet(name="Affymetrix_2006-CEPH270", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2005-500k_data", minCount=5);

if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(dataPath);
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalize
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("qn", mode="list")) {
  verbose && enter(verbose, "Calibration all data sets");
  tmp <- raw;
  for (name in names(raw)) {
    verbose && enter(verbose, "Data set: ", name);
    ds <- raw[[name]];
    path <- findCelSet(name=getName(ds), chipType=getChipType(getCdf(ds)), paths="calibAllelicCT", minCount=length(ds));
    if (is.null(path)) {
      dsN <- calibrateAllelicCrosstalk(ds, verbose=verbose);
    } else {
      verbose && enter(verbose, "Detected calibrated dataset");
      verbose && cat(verbose, "Pathname: ", path);
      dsN <- fromFiles(ds, path);
      verbose && exit(verbose);
    }
    print(dsN);
    tmp[[name]] <- dsN;
    verbose && exit(verbose);
  }
  qn <- tmp;
  verbose && exit(verbose);
}
