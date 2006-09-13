source("init.R");

fig <- 1;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "chip_data3/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup the GDAS annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("gdas")) {
  gdas <- GdasAnnotationSet$forDataSet(ds);
  loci <- join(gdas, orderBy=c("Chromosome", "PhysicalPosition"));
}
chromosomes <- c(1:22, "X");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a set of models to work with
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
print(model);

# Always re-grab the data set, chip effects, and probe affinities.
ds <- getDataSet(model);
samples <- seq(ds);
names <- getNames(ds);

ces <- getChipEffects(model);
print(ces);
cesAvg <- getAverageFile(ces, verbose=TRUE);

ccs <- 23;
samples <- 6;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Display CN estimates along the chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fig <- mm;
fig <- fig + 1;
Device$set(fig, width=12, height=3)
par(ask=FALSE);
for (ss in samples) {
  for (cc in chromosomes[ccs]) {
    cat(sprintf("Sample #%d (%s), chromosome %s...\n", ss, names[ss], cc));
    opar <- par(mar=c(5,4,2,2)+0.1);
    plotMvsPosition(model, sample=ss, chromosome=cc, gdas=gdas, pch=176);
    abline(h=log(1:6/2, base=2), lty=2);
    abline(h=c(-1,1)/2, lty=4);
    par(opar);
    par(ask=FALSE);
  }
}
