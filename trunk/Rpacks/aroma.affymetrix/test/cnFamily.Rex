library(aroma.affymetrix)

pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("dataPath"))
  dataPath <- findCelSet(name="Affymetrix_2006-CEPH270", chipType="Nsp");
print(dataPath)

if (!exists("rawData", mode="list")) {
  verbose && enter(verbose, "Defining all data sets");
  ds <- AffymetrixCelSet$fromFiles(dataPath);
  print(ds);
  raw <- getSiblings(ds);
  rm(ds);
  verbose && exit(verbose);
}
print(raw);


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
  verbose && enter(verbose, "Retrieving reference estimates");
  ceRef <- getAverageFile(ces);
  ceRef$mergeStrands <- model$mergeStrands;
  ceRef$combineAlleles <- model$combineAlleles;
  yRef <- unlist(ceRef[units], use.names=FALSE);
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

  data[[kk]] <- list(x=x, y=y, yRef=yRef, cn=cn, nbrOfSnps=nunits);

  verbose && exit(verbose);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPath <- filePath(getParent(getPath(models[[1]])), "figures");
mkdirs(figPath);

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

for (zoom in c(1,2,4,8,16,32,64,128,256,512)) {
  verbose && enter(verbose, "Plotting zoom x", zoom);

  chrIdx <- which(chr == allChromosomes);
  imgName <- sprintf("%s-MvsX-Chr%02d-x%03d.png", name, chrIdx, zoom);
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
  str <- sprintf("%s=%d SNPs", str, sum(nbrOfSnps));
  stext(str, side=3, pos=1, line=-1, cex=0.8);
  
  # Plot CNs for all chip types
  for (kk in seq(along=data)) {
    df <- data[[kk]];
    x <- df$x;
    cn <- df$cn;
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

