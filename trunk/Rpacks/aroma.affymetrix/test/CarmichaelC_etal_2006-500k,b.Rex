library(aroma.affymetrix)

pngDev <- System$findGraphicsDevice();

chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

dataPath <- "normQuantile";

if (!exists("models", mode="list"))
  models <- list();

for (chipType in chipTypes) {
  model <- models[[chipType]];
  if (!is.null(model))
    next;

  verbose && enter(verbose, "Define combined dataset");

  # Test set
  dataSetName <- "CarmichaelC_etal_2006-500k,b";
  path <- filePath(dataPath, dataSetName, chipType);
  print(path);
  ds <- AffymetrixCelSet$fromFiles(path);
  print(ds);

  # Append reference set
  dataSetName <- "Affymetrix_2006-HapMap270,b";
  path <- filePath(dataPath, dataSetName, chipType);
  print(path);
  dsR <- AffymetrixCelSet$fromFiles(path);
  print(dsR);

  # Append reference set
  append(ds, dsR);
  verbose && exit(verbose);

  # Define model
  verbose && enter(verbose, "Define model");
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
  rm(path, dataSetName, ds, dsR, model, ces, paf);
} # for (chipType ...)
gc()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit copy number model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
allChromosomes <- c(1:22, "X");
#allChromosomes <- c(11, 14, 21);
fitModels <- TRUE;
#fitModels <- FALSE;
if (fitModels) {
  for (chr in allChromosomes) {
    verbose && enter(verbose, "Chromosome ", chr);
  
  for (kk in seq(along=models)) {
    model <- models[[kk]];
    ds <- getDataSet(model);
    verbose && enter(verbose, "Data set: ", getName(ds));
  
    verbose && enter(verbose, "Retrieving genome information");
    cdf <- getCdf(model);
    gi <- getGenomeInformation(cdf);
    units <- getUnitIndices(gi, chromosome=chr);
    nunits <- length(units);
    verbose && exit(verbose);
  
    verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
    uDone <- fit(model, units=units, moreUnits=30, verbose=verbose);
    nDone <- length(uDone);
  
    # Update average chip effects
    if (nDone > 0) {
      ces <- getChipEffects(model, verbose=verbose);
      ceAvg <- getAverageFile(ces, units=uDone, moreCells=30, force=TRUE, verbose=less(verbose));
      ceAvg$mergeStrands <- model$mergeStrands;
      ceAvg$combineAlleles <- model$combineAlleles;
    }
  
    verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
    verbose && exit(verbose);
    verbose && exit(verbose);
  }
    gc();
    verbose && exit(verbose);
  } # for (chr ...)
} # if (fitModels)


stop()

for (chr in allChromosomes) {
  verbose && enter(verbose, "Chromosome ", chr);

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

#sample <- "1-288 [a-zA-Z]* 500K";
sample <- "2-437 [a-zA-Z]* 500K";
#sample <- "3-575 [a-zA-Z]* 500K";

data <- vector("list", length(models));
names(data) <- names(models);
for (kk in seq(along=models)) {
  model <- models[[kk]];
  ds <- getDataSet(model);
  verbose && enter(verbose, "Data set: ", getName(ds));
  ces <- getChipEffects(model);

  verbose && enter(verbose, "Retrieving genome information");
  cdf <- getCdf(model);
  gi <- getGenomeInformation(cdf);
  units <- getUnitIndices(gi, chromosome=chr);
  nunits <- length(units);
  # Get the positions of all SNPs
  x <- getPositions(gi, units=units);
  stopifnot(identical(length(x), nunits));
  verbose && exit(verbose);

  # Remove NAs (SNPs for which we have no position information)
  keep <- !is.na(x);
  x <- x[keep];
  units <- units[keep];
  nunits <- length(units);

  # Get reference (average chip effect)
  verbose && enter(verbose, "Retrieving average reference estimates");
  ceRef <- getAverageFile(ces);
#  ceRef <- getFile(ces, 3);
  ceRef$mergeStrands <- model$mergeStrands;
  ceRef$combineAlleles <- model$combineAlleles;
  yRef <- unlist(ceRef[units], use.names=FALSE);
  sdRef <- readUnits(ceRef, units=units, readIntensities=FALSE, readStdvs=TRUE)
  sdRef <- unlist(sdRef, use.names=FALSE);
  stopifnot(identical(length(yRef), nunits));
  verbose && exit(verbose);

  # Get sample of interest
  verbose && enter(verbose, "Retrieving sample estimates");
  ss <- grep(sample, getNames(ces));
  if (length(ss) == 0)
    throw("No such sample in this dataset: ", sample);
  ce <- getFile(ces, ss);
  verbose && cat(verbose, "Sample:");
  verbose && print(verbose, ce);


  y <- unlist(ce[units], use.names=FALSE);
  stopifnot(identical(length(y), nunits));
  verbose && exit(verbose);

  cn <- y/yRef;
  stopifnot(identical(length(cn), nunits));

  data[[kk]] <- list(x=x, y=y, yRef=yRef, sdRef=sdRef, cn=cn, nbrOfSnps=nunits);

  verbose && exit(verbose);
}

stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPath <- filePath(getParent(getPath(models[[1]])), "figures");
mkdirs(figPath);

# Get genome data
genome <- readTable("annotations/hgChromosomes.txt", colClasses=c(nbrOfBases="integer", nbrOfGenes="integer"), header=TRUE, row.names=1);


name <- grep(sample, getNames(ds), value=TRUE);
nbrOfBases <- genome[chr,"nbrOfBases"];
chipTypes <- sapply(qn, FUN=function(ds) getChipType(getCdf(ds)));
nbrOfSnps <- sapply(data, FUN=function(df) df$nbrOfSnps);

xScale <- 1/10^6;
xlab <- "Physical position (Mb)"
ylab <- expression(log[2](theta/theta[avg]));
xlim <- c(0,nbrOfBases)*xScale;
ylim <- c(-1,1)*2;

pch <- 19;
cex <- 0.8;
cols <- brewer.pal(8, "Dark2");
cols <- rep(cols, length.out=length(data));

zooms <- c(1,2,4,8,16,32,64,128,256,512);
for (zoom in zooms) {
  verbose && enter(verbose, "Plotting zoom x", zoom);

  chrIdx <- which(chr == allChromosomes);
  imgName <- sprintf("%s_b-MvsX-Chr%02d-x%03d.png", name, chrIdx, zoom);
  imgName <- gsub(" ", "_", imgName);
  pathname <- filePath(figPath, imgName);
  width <- round(zoom*5*(nbrOfBases / 10^6));
  height <- 400;
  verbose && printf(verbose, "Image: %s (%.0f,%0.f)\n", 
                                            pathname, width, height);
  pngDev(pathname, width=width, height=400);
  
  opar <- par("mar");
  par(mar=c(5,4,2,2)+0.1);
  
  # Create plot
  plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
  
  # Annotate plot
  abline(h=0, col="blue", lwd=2);
  abline(h=log(1:4/2, base=2), col="blue");
  abline(h=c(-1,1)/2, col="blue", lty=2);
  
  str <- paste(chipTypes, collapse=" / ");
  stext(str, side=4, pos=1, line=0, cex=0.8, col="darkgray");
  stext(sprintf("Sample: %s", sample), side=3, pos=0, line=0, cex=0.8);
  stext(sprintf("Chromosome %s", chr), side=3, pos=1, line=0, cex=0.8);
  str <- paste(nbrOfSnps, collapse="+");
  str <- sprintf("%s=%d SNPs", str, sum(unlist(nbrOfSnps)));
  stext(str, side=3, pos=1, line=-1, cex=0.8);
  
  # Plot CNs for all chip types
  for (kk in seq(along=data)) {
    df <- data[[kk]];
    x <- df$x;
    cn <- df$cn;
    qc <- 1/log2(df$sdRef)^3;
    qc <- qc / median(qc);
#    cex <- 0.2*3/2*qc;
    points(x*xScale, log(cn,base=2), pch=pch, cex=cex, col=cols[kk]);
  }
  rm(df, x, cn);

  par(opar);
  
  # Close and flush figure
  dev.off();

  verbose && exit(verbose);
} # for (zoom in ...)

rm(data);
gc();

  verbose && exit(verbose);
} # for (chr in c(1:22, "X"))

