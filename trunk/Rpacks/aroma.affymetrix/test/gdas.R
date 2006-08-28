source("init.R");

fig <- 1;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup the GDAS annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("gdas")) {
  gdas <- GdasAnnotationSet$fromFiles(path="annotations");
  join(gdas, orderBy=c("Chromosome", "PhysicalPosition"));
}
chromosomes <- c(1:22, "X");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "chip_data3/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a set of models to work with
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  models <- list(
    liwong = AffymetrixTotalCnLiWongModel(ds),
    rma    = AffymetrixTotalCnRmaModel(ds),
    affine = AffymetrixTotalCnAffineModel(ds)
  )
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the model we want to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
mm <- 2;
model <- models[[mm]];
print(model);

# Always re-grab the data set, chip effects, and probe affinities.
ds <- getDataSet(model);
ces <- getChipEffects(model);
samples <- 1:length(ces);

print(ces);

cesAvg <- getAverageFile(ces, verbose=TRUE);

ccs <- 5;
samples <- 19;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Display CN estimates along the chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fig <- mm;
fig <- fig + 1;
Device$set(fig, width=12, height=3)
par(ask=TRUE);
for (ss in samples) {
  ce <- getFile(ces, ss);
  print(getName(ce));
  for (cc in chromosomes[ccs]) {
    cat(sprintf("Sample %d, chromosome %s...\n", ss, cc));
    opar <- par(mar=c(5,4,2,2)+0.1);
    plotMvsPosition(ce, cesAvg, gdas=gdas, chromosome=cc, pch=176);
    abline(h=log(1:6/2, base=2), lty=2);
    abline(h=c(-1,1)/2, lty=4);
    stext(side=3, pos=1, sprintf("Chromosome %s", cc));
    stext(side=1, pos=1, line=-1, cex=0.7, getLabel(model), col="darkgrey");
    par(opar);
  }
}
