source("init.R");

# Specify the dataset to be used
path <- "chip_data/Hind/";
path <- "chip_data2/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# Specify the model we want to fit
model <- AffymetrixCnRmaModel(ds, name="modelTCNRMA");
print(model);

ces <- getChipEffects(model);
cesAvg <- getAverageFile(ces, verbose=TRUE);

if (!exists("gdas")) {
  gdas <- GdasAnnotationSet$fromFiles(path="annotations");
  join(gdas, orderBy=c("Chromosome", "PhysicalPosition"));
}


fig <- 1;

ccs <- c(1:22,"X");
ccs <- c("X");
Device$set(2, width=12, height=3)
par(ask=TRUE);
for (ss in 6:length(ces)) {
  ce <- getFile(ces, ss);
  for (cc in ccs) {
    cat(sprintf("Sample %d, chromosome %s...\n", ss, cc));
    plotMvsPosition(ce, cesAvg, gdas=gdas, chromosome=cc, pch=176);
    abline(h=log(1:6/2, base=2), lty=3);
    stext(side=3, pos=1, sprintf("Chromosome %s", cc));
  }
}
