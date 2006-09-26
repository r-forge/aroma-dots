pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("dataPath"))
  dataPath <- findCelSet(chipType="Xba", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2006-CEPH270", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2005-500k_data", minCount=5);
dataPath <- findCelSet(name="Affymetrix_2004-100k_trios", minCount=5);
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
#    model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
    model <- RmaCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
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

chr <- "1";

allChromosomes <- c(1:22, "X");
chrIdx <- which(chr == allChromosomes);
imgName <- sprintf("%s-SdVsFragLength_Chromosome%02d-2400bp.png", getName(ds), chrIdx);
pathname <- filePath(figPath, imgName);
height <- 640;
width <- 1.618*height / 1.08;
pngDev(pathname, width=width, height=height);

NORMAL <- 1; BOLD <- 2; ITALIC <- 3; BOLDITALIC <- 4;
par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,4,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);

for (kk in seq(along=models)) {
  model <- models[[kk]];
  ds <- getDataSet(model);
  verbose && enter(verbose, "Data set: ", getName(ds));
  ces <- getChipEffects(model);

  verbose && enter(verbose, "Retrieving SNP information");
  cdf <- getCdf(model);
  gi <- getGenomeInformation(cdf);
  si <- getSnpInformation(cdf);
  units <- getUnitIndices(gi, chromosome=chr);
  nunits <- length(units);
  verbose && printf(verbose, "Identified %d SNPs on chromosome %s\n", nunits, chr);
 
  fl <- getFragmentLengths(si, units=units);
  keep <- which(!is.na(fl));
  fl <- fl[keep];
  units <- units[keep];
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
  yRef <- readUnits(ceRef, units=units, verbose=verbose);
  yRef <- unlist(yRef, use.names=FALSE);
  stopifnot(identical(length(yRef), nunits));
  verbose && exit(verbose);

  # Get sample (average chip effect)
  verbose && enter(verbose, "Retrieving sample estimates");
  ss <- grep(sample, getNames(ces));
  if (length(ss) == 0)
    throw("No such sample in this dataset: ", sample);
  ce <- getFile(ces, ss);
  verbose && cat(verbose, "Sample:");
  verbose && print(verbose, ce);
  yS <- readUnits(ce, units=units, verbose=verbose);
  yS <- unlist(yS, use.names=FALSE);
  stopifnot(identical(length(yS), nunits));
  verbose && exit(verbose);

  flMax <- ceiling(max(fl)/100)*100;
  flMin <- floor(min(fl)/100)*100;
  cutOffs <- as.integer(seq(from=flMax, to=flMin, by=-10));
  nCutOffs <- length(cutOffs);
  stats <- matrix(NA, nrow=nCutOffs, ncol=4);
  colnames(stats) <- c("cutoff", "n", "mad", "sd");
  
  bw <- as.integer(300);
  for (cc in seq(along=cutOffs)) {
    cutOff <- cutOffs[cc];
    rr <- cutOff + c(-1,1)*bw/2;
    verbose && enter(verbose, sprintf("Filtering to get SNPs with fragment lengths [%d,%d]", rr[1], rr[2]));
    keep <- which((rr[1] <= fl) & (fl <= rr[2]));
    nkeep <- length(keep);
    verbose && printf(verbose, "Found %d such SNPs\n", nkeep);
    if (nkeep == 0)
      break;
    ySf <- yS[keep]; 
    yRf <- yRef[keep]; 
    cn <- ySf / yRf;
    lcn <- log(cn, base=2);
    mad <- mad(lcn, na.rm=TRUE);
    sd <- sd(lcn, na.rm=TRUE);
    verbose && printf(verbose, "MAD: %.3f, SD: %.3f", mad, sd);
    stats[cc,"cutoff"] <- cutOff;
    stats[cc,"n"] <- nkeep;
    stats[cc,"mad"] <- mad;
    stats[cc,"sd"] <- sd;
    verbose && exit(verbose);
  } # for (cc in ...)
  print(stats);

  xlab <- sprintf("PCR fragment length (window size %d bases)", bw);
  ylab <- expression(hat(sigma));
  xlim <- c(0, 2400);

  ylim <- range(stats[,"mad"], na.rm=TRUE);
  ylim[1] <- 0.10; 
  ylim[2] <- 0.22; 

  if (kk == 1) {
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
    # Annotate plot
#    stext(getName(ces), side=3, pos=0, cex=0.8);
    stext(sprintf("Sample: %s", sample), side=3, pos=0, line=0, cex=0.8);
    stext(sprintf("Chromosome %s", chr), side=3, pos=1, line=0, cex=0.8);
  }
  keep <- (stats[,"n"] > 15);
  lines(stats[keep,c("cutoff", "mad")], col=kk, lwd=10);

  # Add density of fragment lengths
  usr <- par("usr");
  d <- density(fl, adjust=1/5);
  d$y <- usr[3] + d$y / max(d$y) * (diff(usr[3:4]) / 10);
  lines(d, col=kk, lwd=2);
} # for (kk in seq(along=models))

ct <- sapply(models, FUN=function(m) getChipType(getCdf(m)));
legend(x="topleft", legend=ct, col=1:2, pch=19, cex=1, bty="n", y.intersp=1.5);

dev.off();
