library(aroma.affymetrix);
source("init.R");


plotSdEsts <- function(all, xlim=c(1,42), ylim=c(0,0.15), ...) {
  opar <- par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,2,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);
  on.exit(par(opar));
  
  plot(NA, xlim=xlim, ylim=ylim, xlab="#samples in reference set", ylab=expression(hat(sigma)));

  for (ss in 1:nrow(all[[1]])) {
    lC <- lapply(all, FUN=function(b) b[ss,]);
    lC <- as.data.frame(lC);
    lC <- as.matrix(lC);
    lC <- unname(lC);

    sdEst <- apply(lC, MARGIN=1, FUN=mad);
    x <- 1:length(sdEst);
#    points(x, sdEst, cex=1, lwd=1, col="lightgray");
    fit <- robustSmoothSpline(x, sdEst, spar=0.8);
    lines(fit, col=ss, lwd=3)
  }
}

#stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pngDev <- System$findGraphicsDevice();
ext <- "png";
if (ext == "png") {
  height <- 640;
  dev <- pngDev;
} else if (ext == "eps") {
  height <- 8 / 1.618;
  dev <- eps;
}
#dev <- NULL;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The test set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("ces")) {
  path <- "modelRmaCnPlm/Affymetrix_2006-HapMap270,HapMap,QN,CEU,RMA,A+B/Mapping250K_Nsp";
#  path <- "modelRmaCnPlm/CarmichaelC_etal_2006-500k,AGRFall,time,nref35,QN,RMA,A+B/Mapping250K_Nsp";
#  path <- "normFragmentLength/Affymetrix_2006-HapMap270,HapMap,QN,CEU,RMA,A+B,FLN/Mapping250K_Nsp";
#  path <- "normFragmentLength/CarmichaelC_etal_2006-500k,AGRFall,time,nref35,QN,RMA,A+B,FLN/Mapping250K_Nsp";
  ces <- CnChipEffectSet$fromFiles(path);
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);
  rm(anno);
}

if (!exists("anno")) {
  cdf <- getCdf(ces);
  gi <- getGenomeInformation(cdf);
  uu <- getUnitsOnChromosome(gi, 3);
  pos <- getPositions(gi, units=uu);
#  library(gtools); # running()
#  dw <- running(pos, width=5, fun=function(x) diff(range(x)));
  
  units <- c(255379, 63612, 63608, 250296, 243393);

#  units <- c(234891, 33613, 221345, 20044, 238717);
  unitNames <- getUnitNames(cdf, units=units);
  anno <- getData(gi, units=units);
  anno <- cbind(snp=unitNames, anno);
  if (length(unique(anno$chromosome)) != 1)
    throw("Internal error");
  print(anno);
  rm(yUnits);
}

# Get data
if (!exists("yUnits")) {
  yUnits <- readUnits(ces, units=units);
}

# Extract data as a matrix
y <- lapply(yUnits, FUN=function(unit) as.numeric(unit[[1]]$theta));
y <- as.data.frame(y);
rownames(y) <- getNames(ces);
colnames(y) <- gsub("[.]", "-", colnames(y));
y <- t(y);

snps <- seq(nrow(y));
#snps <- 1;
y <- y[snps,,drop=FALSE];

  ylim <- c(-0.3,0.7);
  ylim <- c(-0.4,0.2);

# Extract one sample
ss <- which(getNames(ces) == "1-288");
sampleName <- colnames(y)[ss];
yS <- y[,ss];
lyS <- log2(yS);

# Use rest as a reference
#y <- y[,-ss,drop=FALSE];
#keep <- !(getNames(ces) %in% c("1-288", "2-437", "3-575")); y <- y[,keep,drop=FALSE];
ly <- log2(y);

chipType <- getChipType(cdf);
chipType <- gsub("-monocell$", "", chipType);


Bs <- c(1,2,10,30,200);
for (B in as.integer(Bs)) {
  set.seed(1);

  bbs <- 1:B;

  tags <- NULL;
  tags <- c(tags, setdiff(getTags(ces), c("HapMap", "QN", "A+B")));
  tags <- c(tags, sprintf("%dSNPs", length(yS)), "MvsRefSize");
  tags <- c(tags, sprintf("%04drun",B));
  if (B <= 3) {
    lwd=3;
  } else {
    lwd=1;
  }
  pathname <- sprintf("%s,%s.%s", sampleName, paste(tags, collapse=","), ext);
  width <- height / 1.08;
  if (!is.null(dev)) {
    print(pathname);
    dev(pathname, width=width, height=height);
  }
  
  opar <- par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,2,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);
  
  xlim <- c(1,ncol(y));
  xlim <- c(1,45);
  plot(NA, xlim=xlim, ylim=ylim, xlab="#samples in reference set", ylab=expression(log[2](theta)-log[2](theta[R])));
  abline(h=0, lty=2, lwd=2, col="lightgray");
  stext(sprintf("Sample: %s", sampleName), side=3, pos=0, line=0.2, cex=0.8);
  if (length(snps) > 1)
    stext(sprintf("SNP range: %d bases", sum(diff(anno$physicalPosition[snps]))), side=3, pos=1, line=0.2, cex=0.8);
  stext(chipType, side=4, pos=1, line=0.4, cex=0.7, col="gray");
  if (B > 1)
    stext(sprintf("#runs: %d", B), side=3, pos=1, line=-1.4, cex=0.7);
  
  all <- list();
  for (bb in bbs) {
    lyAvg <- matrix(NA, nrow=nrow(y), ncol=ncol(y));
    fullSubset <- sample(ncol(y));
    for (kk in seq(from=1, to=ncol(y))) {
      subset <- fullSubset[1:kk];
      subset <- sample(ncol(y), size=kk);
      lyAvg[,kk] <- apply(ly[,subset,drop=FALSE], MARGIN=1, FUN=median, na.rm=TRUE);
    }
    lC <- lyS - lyAvg; 
    all[[bb]] <- lC;
  }

  for (ss in 1:nrow(y)) {
    lC <- lapply(all, FUN=function(b) b[ss,]);
    lC <- as.data.frame(lC);
    lC <- as.matrix(lC);
    lC <- unname(lC);

    apply(lC, MARGIN=2, FUN=lines, lwd=lwd, col=ss);
  }

  par(opar);
  if (!is.null(dev))
    dev.off();
} # for (B ...)


# Plot SD estimate from last run
tags <- c(tags[-length(tags)], "SDvsRefSize");
pathname <- sprintf("%s,%s.%s", sampleName, paste(tags, collapse=","), ext);
width <- height / 1.08;
if (!is.null(dev)) {
  print(pathname);
  dev(pathname, width=width, height=height);
}
opar <- par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,2,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);
 
plotSdEsts(all, xlim=xlim, ylim=c(0,0.16));
stext(sprintf("Sample: %s", sampleName), side=3, pos=0, line=0.2, cex=0.8);
stext(chipType, side=4, pos=1, line=0.4, cex=0.7, col="gray");
if (length(bbs) > 1)
  stext(sprintf("#runs: %d", B), side=3, pos=1, line=-1.4, cex=0.7);
par(opar);
if (!is.null(dev))
  dev.off();

