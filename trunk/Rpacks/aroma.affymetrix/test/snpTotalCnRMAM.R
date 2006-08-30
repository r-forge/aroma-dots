source("init.R")

verbose <- Arguments$getVerbose(TRUE);
fig <- 1;
png <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
#path <- "chip_data/Hind/";
path <- "SlaterH_2004/Xba/";
#path <- "chip_data2/Sty/";

if (!exists("ds")) {
  ds <- AffymetrixCelSet$fromFiles(path);
  chipType <- getChipType(getCdf(ds))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup the GDAS annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("gdas")) {
#  gdas <- GdasAnnotationSet$fromFiles(path=filePath("annotations", chipType));
#  join(gdas, orderBy=c("Chromosome", "PhysicalPosition"));
}
allChromosomes <- c(1:22, "X");
chromosomes <- allChromosomes;


if (!exists("model")) {
  # Specify the model we want to fit
  model <- AffymetrixTotalCnRmaModel(ds);
  print(model);
  rm(ds);  # Not need anymore
}

# Always re-grab the data set, chip effects, and probe affinities,
# after setting up a model.
ds <- getDataSet(model);
samples <- seq(ds);
names <- getNames(ds);

paf <- getProbeAffinities(model);
ces <- getChipEffects(model);
print(paf);
print(ces);

# Test to fit a single unit first
#units <- fit(model, units=56, force=TRUE, verbose=TRUE);

#units <- findUnitsTodo(paf);
#print(summary(units));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- fit(model, moreUnits=10, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculate the average copy-number file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cesAvg <- getAverageFile(ces, verbose=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot M vs A
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)
  ce <- getFile(ces, 1);
  smoothScatterMvsA(ce, cesAvg, xlim=c(8,16));
  abline(h=log(1:6/2, base=2), lty=3);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Display CN estimates along the chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chromosomes <- "5";
samples <- 19;

cdf <- getCdf(model);
figPath <- filePath("figures", getChipType(cdf));
mkdirs(figPath);

offset <- -300;

fig <- fig + 1;
for (ss in samples) {
  ce <- getFile(ces, ss);
  ce$offset <- offset; 

  name <- getLabel(ce);
  imgname <- sprintf("%s-MvsA-offset%+03d.png", name, offset);
  pathname <- filePath(figPath, imgname);  
  if (!isFile(pathname)) {
    png(pathname, width=800, height=800);
    smoothScatterMvsA(ce, cesAvg, xlim=c(8,16));
    abline(h=log(1:6/2, base=2), lty=3);
    dev.off();
  }
  
  for (cc in chromosomes) {
    cat(sprintf("Sample #%d (%s), chromosome %s...\n", ss, names[ss], cc));

    tt <- select(gdas, Chromosome=cc, sortBy="integer:PhysicalPosition");
    tt <- as.integer(tt[c(1,nrow(tt)),"PhysicalPosition"]);
    chrLength <- tt[2]/1e6;
    
    imgname <- sprintf("%s-Chr%02d-CN-offset%+03d.png", name, which(cc == allChromosomes), offset);
    pathname <- filePath(figPath, imgname);  
    if (!isFile(pathname)) {
      png(pathname, width=20*chrLength, height=400);
      opar <- par(mar=c(5,4,2,2)+0.1);
      plotMvsPosition(model, sample=ss, chromosome=cc, gdas=gdas, pch=19, ylim=c(-1,1)*1.5);
      abline(h=log(1:6/2, base=2), lty=2);
      par(opar);
      dev.off();
    }

    imgname <- sprintf("%s-Chr%02d-A-offset%+03d.png", name, which(cc == allChromosomes), offset);
    pathname <- filePath(figPath, imgname);  
    if (!isFile(pathname)) {
      png(pathname, width=20*chrLength, height=400);
      opar <- par(mar=c(5,4,2,2)+0.1);
      plotMvsPosition(model, sample=ss, chromosome=cc, gdas=gdas, pch=19, ylim=c(8,16), what="A");
      abline(h=seq(from=0, to=16, by=2), lty=2);
      par(opar);
      dev.off();
    }
  }
}

