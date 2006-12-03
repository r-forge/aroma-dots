library(aroma.affymetrix);
source("init.R");

ans <- as.integer(readline("Enter 1-4: "));
print(ans);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("dsHapMap")) {
  rootPath <- filePath("raw/Affymetrix_2006-HapMap270/");
  pathname <- filePath(rootPath, "HapMap,origin.txt");
  samples <- read.table(file=pathname, colClasses=rep("character", 2), header=TRUE, sep="\t");

  rootPath <- filePath("normQuantile/Affymetrix_2006-HapMap270,QN/");
  path <- filePath(rootPath, "Mapping250K_Nsp");
  ds <- AffymetrixCelSet$fromFiles(path);

  # Reorder samples by according HapMap list/groups.
  reorder(ds, order=samples$sample);
  date <- format(getTimestamps(ds), "%Y-%m-%d");
  samples <- cbind(samples, date=I(date));

  dsT <- list();

  # The YRI samples where hybridized mostly on 2005-12-06
  keep <- (samples$origin == "YRI") & (date == "2005-12-06");
  keep <- which(keep)[1:40];
  yri <- extract(ds, keep);
  setTags(yri, c(getTags(yri), "YRI"));
  dsT$YRI <- yri;

  # The CHB samples where hybridized mostly on 2005-11-11 & 2005-11-22
  keep <- (samples$origin == "CHB") & (date == "2005-11-11" | date == "2005-11-22");
  keep <- which(keep)[1:40];
  chb <- extract(ds, keep);
  setTags(chb, c(getTags(chb), "CHB"));
  dsT$CHB <- chb;

  # The CEU samples where hybridized mostly on 2005-11-18 & 2005-11-19
  keep <- (samples$origin == "CEU") & (date == "2005-11-18" | date == "2005-11-19");
  keep <- which(keep)[1:40];
  ceu <- extract(ds, keep);
  setTags(ceu, c(getTags(ceu), "CEU"));
  dsT$CEU <- ceu;

  # The JPT samples where hybridized mostly on 2005-11-22 & 2005-11-23
  keep <- (samples$origin == "JPT") & (date == "2005-11-22" | date == "2005-11-23");
  keep <- which(keep)[1:40];
  jpt <- extract(ds, keep);
  setTags(jpt, c(getTags(jpt), "JPT"));
  dsT$JPT <- jpt;

  dsHapMap <- dsT;
  rm(ds, dsT);
  gc();
}

print(dsHapMap);

if (!exists("dsT")) {
  path <- "normQuantile/CarmichaelC_etal_2006-500k,HapMap,QN/Mapping250K_Nsp";
  dsT <- AffymetrixCelSet$fromFiles(path);
}
print(dsT);


res <- list();
for (kk in seq(along=dsHapMap)[ans]) {
  ds <- clone(dsT);
  append(ds, dsHapMap[[kk]]);
  setName(ds, getName(dsHapMap[[kk]]));
  setTags(ds, unique(c(getTags(ds), getTags(dsHapMap[[kk]]))));
  print(ds);  
  gc();

  # Fit PLM
  verbose && enter(verbose, "Summarazing probe signals");
  plm <- RmaCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);

  verbose && print(verbose, plm);
  fit(plm, moreUnits=1, verbose=verbose);
  ces <- getChipEffects(plm);
  rm(plm);
  gc();
  verbose && exit(verbose);

  # Normalize for PCR fragment-length effects
  verbose && enter(verbose, "Normalizing summarized SNP signals");
  verbose && print(verbose, ces);
  normFL <- FragmentLengthNormalization(ces, subsetToFit=1/3);
  rm(ces);
  verbose && print(verbose, normFL);
  cesFL <- process(normFL, verbose=verbose);
  rm(normFL);
  gc();
  verbose && print(verbose, cesFL);
  verbose && exit(verbose);

  gc();
}
