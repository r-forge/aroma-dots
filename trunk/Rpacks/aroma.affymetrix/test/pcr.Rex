library(aroma.affymetrix);
source("init.R");

pngDev <- System$findGraphicsDevice();
ext <- "png";
if (ext == "png") {
  height <- 640;
  dev <- pngDev;
} else if (ext == "eps") {
  height <- 8 / 1.618;
  dev <- eps;
}

chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

dataPath <- "normQuantile";

if (!exists("models", mode="list"))
  models <- list();

for (chipType in chipTypes) {
  model <- models[[chipType]];
  if (!is.null(model))
    next;

  verbose && enter(verbose, "Define dataset");

  # Test set
  dataSetName <- "SlaterH_etal_2004,b";
  path <- filePath(dataPath, dataSetName, chipType);
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  verbose && exit(verbose);

  # Define model
  verbose && enter(verbose, "Define model");
#  model <- RmaCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
  model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);

  verbose && enter(verbose, "Pre-creates probe affinities");
  paf <- getProbeAffinities(model);
  verbose && exit(verbose);

  verbose && enter(verbose, "Pre-creates chip effects");
  ces <- getChipEffects(model);
  verbose && exit(verbose);

  verbose && exit(verbose);

  models[[chipType]] <- model;
 
  # Not needed anymore
  rm(path, dataSetName, ds, model, ces, paf);
} # for (chipType ...)
gc()

doRatio <- F;
par(ask=FALSE);
chr <- 2;
chr <- NULL;
chrIdx <- chr;

xlab <- "PCR fragment length (bases)";
xlim <- c(0,2400);
if (doRatio) {
  ylab <- expression(log[2](C));
  ylim <- c(-1,1)*4;
} else {
  ylab <- expression(log[2](theta));
  ylim <- c(8,18);
}

for (kk in seq(along=models)) {
  model <- models[[kk]];
  ds <- getDataSet(model);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);
  si <- getSnpInformation(cdf);
  ces <- getChipEffects(model);
  ceAvg <- getAverage(ces);
  ceAvg$combineAlleles <- model$combineAlleles;
  ceAvg$mergeStrands <- model$mergeStrands;

  verbose && enter(verbose, "Data set: ", getName(ces));
  verbose && cat(verbose, "Chip type: ", getChipType(cdf));

  figPath <- filePath("pcrFragmentLength", getFullName(ces), getChipType(cdf));
  mkdirs(figPath);

  verbose && enter(verbose, "Indentifying units");
  if (is.null(chr)) {
    units <- getUnitIndices(gi);
    units <- units[seq(from=1, to=length(units), length=10000)];
  } else {
    units <- getUnitIndices(gi, chromosome=chr);
  }
  fl <- getFragmentLengths(si, units=units);
  keep <- !is.na(fl);
  fl <- fl[keep];
  units <- units[keep];
  verbose && str(verbose, units);
  verbose && exit(verbose);

  # Density of fragment lengths
  d <- density(fl, adjust=1/2);
  d$y <- ylim[1] + d$y / max(d$y) * (diff(ylim) / 10);
  flD <- d;
 
  verbose && enter(verbose, "Retrieving reference signals");
  yR <- readUnits(ceAvg, units=units);
  yR <- unlist(yR, use.names=FALSE);
  verbose && str(verbose, yR);
  verbose && exit(verbose);
  
  for (ss in seq(ces)) {
    ce <- getFile(ces, ss); 
    verbose && enter(verbose, sprintf("Sample #%d (%s)", ss, getName(ce)));

    if (is.null(chr)) {
      imgName <- sprintf("%s", getName(ce));
    } else {
      imgName <- sprintf("%s,Chr%02d", getName(ce), chrIdx);
    }
    if (doRatio) {
      imgName <- sprintf("%s,CnVsFragLength", imgName);
    } else {
      imgName <- sprintf("%s,ThetaVsFragLength", imgName);
    }
    imgName <- sprintf("%s.%s", imgName, ext);
    pathname <- filePath(figPath, imgName);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already plotted...skipped");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Retrieving signals");
    y <- readUnits(ce, units=units);
    y <- unlist(y, use.names=FALSE);
    verbose && str(verbose, yR);
    verbose && exit(verbose);

    if (doRatio) {
      lc <- log2(y/yR);
    } else {
      lc <- log2(y);
    }

    width <- 1.618*height / 1.08;
    dev(pathname, width=width, height=height);

    NORMAL <- 1; BOLD <- 2; ITALIC <- 3; BOLDITALIC <- 4;
    par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,2,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);

    sample <- sprintf("#%d. %s", ss, getName(ce));
    smoothScatter(fl, lc, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=NULL);
    abline(h=0, lwd=4, col="orange");
    ok <- is.finite(lc);
    l <- lowess(fl[ok], lc[ok]);
    lines(l, lwd=4, col="red");
    lines(flD, col="darkgray", lwd=2);

    # Annotate plot
    stext(sprintf("Sample: %s", sample), side=3, pos=0, line=0, cex=0.8);
    stext(sprintf("%d SNPs", length(y)), side=3, pos=1, line=-1, cex=0.8);
    stextChipType(cdf, line=0);
    if (!is.null(chr))
      stext(sprintf("Chromosome %s", chr), side=3, pos=1, line=0, cex=0.8);

    dev.off();

    verbose && exit(verbose);
  }

  verbose && exit(verbose);
}
