pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "data";

name <- "Affymetrix_2006-CEPH270";
name <- "Affymetrix_2005-500k_data";
name <- "Affymetrix_2004-100k_trios";
name <- "Affymetrix_2004-100k_trios-calibACT";
name <- "Affymetrix_2005-500k_data-calibACT";
name <- "SinclairA_2006-calibACT";
#name <- "SlaterH_etal_2004-calibACT";

if (!exists("dataPath"))
  dataPath <- findCelSet(name=name, paths=paths, minCount=5);
print(dataPath);

if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(dataPath);
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}
print(raw);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalize
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("qn", mode="list")) {
  verbose && enter(verbose, "Normalizing all data sets");
  tmp <- raw;
  for (name in names(raw)) {
    verbose && enter(verbose, "Data set: ", name);
    ds <- raw[[name]];
    path <- findCelSet(name=getName(ds), chipType=getChipType(getCdf(ds)), paths="normQuantile", minCount=length(ds));
    if (is.null(path)) {
      dsN <- normalizeQuantile(ds, subsetToAvg=1/3, verbose=verbose);
    } else {
      verbose && enter(verbose, "Detected normalized dataset");
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
print(qn);
