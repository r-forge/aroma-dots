pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("dataPath"))
  dataPath <- findCelSet(chipType="Xba", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2006-CEPH270", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2005-500k_data", minCount=5);
#dataPath <- findCelSet(name="Affymetrix_2004-100k_trios", minCount=5);
#dataPath <- findCelSet(name="Affymetrix_2004-100k_trios-calibACT", minCount=5);
#dataPath <- findCelSet(name="Affymetrix_2005-500k_data-calibACT", minCount=5);

if (!exists("raw", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(dataPath);
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}
print(raw);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Normalizing data
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
      print(dsN);
    } else {
      verbose && enter(verbose, "Detected normalized dataset");
      verbose && cat(verbose, "Pathname: ", path);
      dsN <- fromFiles(ds, path);
      verbose && exit(verbose);
    }
    tmp[[name]] <- dsN;
    verbose && exit(verbose);
  }
  qn <- tmp;
  verbose && exit(verbose);
}
print(qn);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create copy-number models
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  verbose && enter(verbose, "Defining PLM models");
  models <- vector("list", length(qn));
  names(models) <- names(qn);
  for (kk in seq(along=models)) {
    ds <- qn[[kk]];
    verbose && enter(verbose, "Data set: ", getName(ds));
    model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
    models[[kk]] <- model;
  }
  verbose && print(verbose, model);

  verbose && enter(verbose, "Pre-creates probe affinities");
  paf <- getProbeAffinities(model, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Pre-creates chip effects");
  ces <- getChipEffects(model, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && exit(verbose);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get relative copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# CEPH 100K family where there is a short amplification in father
# and daughter on Chromosome 8.  See CRLMM poster.
family <- c(mother=NA, father="NA12814", daughter="NA12802");
sample <- family["daughter"];

# The fifthteen CEPH samples that are available in both the 100K set
# (90 samples) and the _public_ 500K set (48 samples)
# > intersect(substr(getNames(ds2), 1,7), substr(getNames(ds), 1,7))
commonCEPH <-c("NA10851", "NA10855", "NA10863", "NA11831", "NA11832", "NA12056", "NA12057", "NA12234", "NA12264", "NA12707", "NA12716", "NA12717", "NA12801", "NA12812", "NA12813")
sample <- commonCEPH[1];


figPath <- filePath(getParent(getPath(models[[1]])), "figures");
mkdirs(figPath);
ext <- "eps";
if (ext == "png") {
  height <- 640;
  dev <- pngDev;
} else if (ext == "eps") {
  height <- 8 / 1.618;
  dev <- eps;
}

for (kk in seq(along=models)) {
  model <- models[[kk]];
  ds <- getDataSet(model);
  verbose && enter(verbose, "Data set: ", getName(ds));
  ces <- getChipEffects(model);
  unitsTodo <- findUnitsTodo(ces, verbose=less(verbose));

  verbose && enter(verbose, "Retrieving SNP information");
  cdf <- getCdf(model);
  si <- getSnpInformation(cdf);
  fl <- getFragmentLengths(si);
  keep <- which(!is.na(fl));
  keep <- setdiff(keep, unitsTodo);
  fl <- fl[keep];
  units <- keep;
  nunits <- length(units);
  stopifnot(identical(length(fl), nunits));
  verbose && printf(verbose, "Identified PCR fragment length details for %d SNPs", nunits);
  verbose && cat(verbose, "Summary:");
  verbose && summary(verbose, fl);
  verbose && exit(verbose);

  # Get reference (average chip effect)
  verbose && enter(verbose, "Retrieving reference estimates");
  ceRef <- getAverageFile(ces);
  ceRef$mergeStrands <- model$mergeStrands;
  ceRef$combineAlleles <- model$combineAlleles;
#  yRef <- unlist(ceRef[units], use.names=FALSE);
  yRef <- readUnits(ceRef, units=units, verbose=verbose);
  yRef <- unlist(yRef, use.names=FALSE);
  stopifnot(identical(length(yRef), nunits));
  verbose && exit(verbose);

  # The fragment-length effect is much more profound on
  # the logartihmic scale.
  log <- TRUE; 
  xlab <- "PCR fragment length (bases)";
  c <- 100;
  xlim <- c(floor(min(fl)/500)*500, ceiling(max(fl)/300)*300);
  xlim <- c(0,2400);
  xlim[1] <- 0;
  y <- yRef;
  if (log) {
    y <- log(y, base=2);
    ylim <- c(6,16);
    ylab <- expression(log[2](theta[ref]));
  } else {
    ylim <- c(0,65535);
    ylab <- expression(theta[ref]);
  }


#  Device$set(kk+1, width="162%", height="108%");
  imgName <- sprintf("%s-%s_ThetaVsFragLength-2400bp.%s", getName(ds), getChipType(cdf), ext);
  pathname <- filePath(figPath, imgName);
  width <- 1.618*height / 1.08;
  dev(pathname, width=width, height=height);
  NORMAL <- 1; BOLD <- 2; ITALIC <- 3; BOLDITALIC <- 4;
  par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,4,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);

  # Smooth scatter plot
  smoothScatter(fl, y, pch=176, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);

  # Annotate plot
  stextChipType(cdf, line=0);
  stext(getName(ces), side=3, pos=0, cex=0.8);
  stext(sprintf("%d samples averaged", length(ces)), side=3, pos=0, line=-1, cex=0.8);
  stext(sprintf("%d SNPs", length(y)), side=3, pos=1, line=-1, cex=0.8);

  # Add density of fragment lengths
  usr <- par("usr");
  d <- density(fl, adjust=1/2);
  d$y <- usr[3] + d$y / max(d$y) * (diff(usr[3:4]) / 10);
  lines(d, col="darkgray", lwd=2);

  # Add lowess line
  ok <- is.finite(y);
  l <- lowess(fl[ok], y[ok]);
  lines(l, col=2, lwd=10);

  dev.off();
}