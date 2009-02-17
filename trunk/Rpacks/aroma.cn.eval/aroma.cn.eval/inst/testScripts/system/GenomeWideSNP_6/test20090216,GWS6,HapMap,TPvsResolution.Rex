if (interactive()) savehistory();
library("aroma.cn.eval");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

tagsList <- list("GTC", "dChip,AVG,ref=median", "dChip,MBEI,ref=median");
dataSet <- "HapMap270,6.0,CEU,founders,three";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load raw CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- lapply(tagsList, FUN=function(tags) {
  AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*");
});
# Keep only common samples (just in case)
names <- Reduce(intersect, lapply(dsList, FUN=getNames));
dsList <- lapply(dsList, FUN=extract, names);
print(dsList);


tags <- Reduce(intersect, lapply(dsList, FUN=getTags));
methods <- sapply(dsList, FUN=function(ds) {
  paste(setdiff(getTags(ds), c(tags, "ref=median")), collapse=",");
})
platforms <- sapply(dsList, FUN=getPlatform);
chipTypes <- sapply(dsList, FUN=getChipType);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Region of interest with known CN states
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample of interest
array <- "NA06985";

# Chromosomal region of interest
chromosome <- 2;
region <- c(81,86)*1e6;

# A function give the true CN state for a given set of positions
truth <- function(x, name=NULL, ...) {
  t <- integer(length(x));
  # Location of change points
  xCps <- c(83.09, 83.72)*1e6;
  # Neutral
  t[x <= xCps[1] | x >= xCps[2]] <- 0L;
  # Loss
  t[xCps[1] <= x & x < xCps[2]] <- -1L;
  # Unknown (safety region of 50kb each side)
  dx <- 50e3;
  for (kk in seq(along=xCps)) {
    xCp <- xCps[kk];
    t[xCp-dx < x & x < xCp+dx] <- NA;
  }
  t;
} # truth()

name <- getNames(extract(dsList[[1]], array))[1];
print(name);
 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract data from region of interest
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cnList <- lapply(dsList, FUN=function(ds) {
  df <- getFile(extract(ds, array), 1);
  cn <- extractRawCopyNumbers(df, chromosome=chromosome, region=region);
  if (is.element("dChip", getTags(ds)))
    cn$cn <- cn$cn-1;
  cn <- SegmentedCopyNumbers(cn=cn$cn, x=cn$x, states=truth);  
  cn;
})
names(cnList) <- sprintf("%s\n%s\n%s", methods, platforms, chipTypes);
print(cnList);
nbrOfTracks <- length(cnList);

# Range of loci covered in all platforms
xRange <- range(sapply(cnList, FUN=xRange));  



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TP rate at fixe FP rate at various resolutions (in kb)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
devSet("TPvResolution");
states <- c(-1,0);
fpRate <- 0.02;
binWidths <- c(0.1,0.5,1,2,4,6,8,10,15,20)*1e3;
xlim <- c(0,max(binWidths))/1e3;
tplim <- c(0,1);
tplab <- "True-positive rate";
xlab <- "Width of bins (in kb)";
par(mar=c(3,3,2,1)+0.1, mgp=c(1.4,0.4,0));
plot(NA, xlim=xlim, ylim=tplim, xlab=xlab, ylab=tplab);
stext(side=3, pos=0, name);
trackCols <- seq(length=nbrOfTracks);
labels <- methods;
legend("bottomright", col=trackCols, pch=19, labels, cex=0.8, bty="n");

naValue <- as.double(NA);
tpRates <- matrix(naValue, nrow=length(binWidths), ncol=nbrOfTracks);
colnames(tpRates) <- names(cnList);
res <- cbind(data.frame(binWidth=binWidths), tpRates);

for (ww in seq(along=binWidths)) {
  binWidth <- binWidths[ww];

  # Smooth CNs using consecutive bins of given width (in kb)
  cnSList <- lapply(cnList, FUN=function(cn) {
    cnS <- binnedSmoothingByState(cn, from=xRange[1], to=xRange[2], by=binWidth, verbose=log);
    cnS <- extractSubsetByState(cnS, states=states);
    cnS;
  })
  print(cnSList);

  fits <- lapply(cnSList, FUN=function(cnS) {
    findRocTpAtFp(cnS, fpRate=fpRate, states=states, recall=states[1]);
  });

  tpRates <- sapply(fits, FUN=function(fit) fit[["tpRateEst"]]);
  cols <- match(names(tpRates), colnames(res));
  row <- whichVector(res[,"binWidth"] == binWidth);
  res[row,cols] <- tpRates;

  x <- rep(binWidth/1e3, nbrOfTracks);
  points(x, y=tpRates, pch=19, col=trackCols);
} # for (ww ...)

o <- order(res[,"binWidth"]);
res <- res[o,,drop=FALSE];

for (kk in seq(length=nbrOfTracks)) {
  binWidths <- res[,"binWidth"];
  tpRates <- res[,kk+1];
  lines(binWidths/1e3, tpRates, lwd=2, col=trackCols[kk]);
}

