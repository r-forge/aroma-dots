savehistory();

figPath <- "figures";
mkdirs(figPath);

pngDev <- System$findGraphicsDevice();

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
# Fit copy number model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get genome data
genome <- readTable("annotations/hgChromosomes.txt", colClasses=c(nbrOfBases="integer", nbrOfGenes="integer"), header=TRUE, row.names=1);


allChromosomes <- c(1:22, "X");
allChromosomes <- allChromosomes[8];
for (chr in allChromosomes) {
  verbose && enter(verbose, "Chromosome ", chr);

if (!exists("models", mode="list")) {
  models <- vector("list", length(qn));
  names(models) <- names(qn);
}

for (kk in seq(along=qn)) {
  ds <- qn[[kk]];
  verbose && enter(verbose, "Data set: ", getName(ds));
  model <- models[[kk]];
  if (is.null(model)) {
    model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
    models[[kk]] <- model;
  }
  verbose && print(verbose, model);

  verbose && enter(verbose, "Retrieving genome information");
  cdf <- getCdf(model);
  gi <- getGenomeInformation(cdf);
  units <- getUnitIndices(gi, chromosome=chr);
  nunits <- length(units);
  verbose && exit(verbose);

  verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", nunits, chr, getChipType(cdf)));
  uDone <- fit(model, units=units, moreUnits=10, verbose=verbose);
  nDone <- length(uDone);

  # Update average chip effects
  if (nDone > 0) {
    ces <- getChipEffects(model);
    ceAvg <- getAverageFile(ces, force=TRUE, verbose=less(verbose));
    ceAvg$mergeStrands <- model$mergeStrands;
    ceAvg$combineAlleles <- model$combineAlleles;
  }

  verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
  verbose && exit(verbose);
  verbose && exit(verbose);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get relative copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sample <- "NA12802";

if (exists("opar"))
  par(opar);

name <- grep(sample, getNames(ds), value=TRUE);
imgName <- sprintf("%s-MvsX.png", name);
pathname <- filePath(figPath, imgName);
nbrOfBases <- genome[chr,"nbrOfBases"];
width <- round(100*(nbrOfBases / 10^6));
pngDev(pathname, width=width, height=400);

opar <- par("mar");
par(mar=c(5,4,2,2)+0.1);
xScale <- 1/10^6;
xlab <- "Physical position (Mb)"
ylab <- expression(log[2](theta/theta[avg]));
xlim <- c(0,nbrOfBases)*xScale;
ylim <- c(-1,1)*2;

pch <- 19;
cex <- 0.8;
cols <- brewer.pal(8, "Dark2");
cols <- rep(cols, length.out=length(qn));

# Create plot
plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main);

# Annotate plot
abline(h=0, col="blue", lwd=2);
abline(h=log(1:4/2, base=2), col="blue");
abline(h=c(-1,1)/2, col="blue", lty=2);

chipTypes <- sapply(qn, FUN=function(ds) getChipType(getCdf(ds)));
str <- paste(chipTypes, collapse=" / ");
stext(str, side=4, pos=1, line=0, cex=0.8, col="darkgray");
stext(sprintf("Sample: %s", sample), side=3, pos=0, line=0, cex=0.8);
stext(sprintf("Chromosome %s", chr), side=3, pos=1, line=0, cex=0.8);

chipTypes <- nbrOfSnps <- c();
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
  ce <- getFile(ces, ss);
  y <- unlist(ce[units], use.names=FALSE);
  stopifnot(identical(length(y), nunits));
  verbose && exit(verbose);

  cn <- y/yRef;
  stopifnot(identical(length(cn), nunits));

  nbrOfSnps <- c(nbrOfSnps, nunits);

  # Plot CNs
  points(x*xScale, log(cn,base=2), pch=pch, cex=cex, col=cols[kk]);

  verbose && exit(verbose);
}
str <- paste(nbrOfSnps, collapse="+");
str <- sprintf("%s=%d SNPs", str, sum(nbrOfSnps));
stext(str, side=3, pos=1, line=-1, cex=0.8);
par(opar);

# Close and flush figure
dev.off();

  verbose && exit(verbose);
} # for (chr in c(1:22, "X"))

