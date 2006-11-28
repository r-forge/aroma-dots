library(aroma.affymetrix);
source("init.R");
source("ui.R");

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
#dev <- NULL; par(ask=TRUE);

NORMAL <- 1; BOLD <- 2; ITALIC <- 3; BOLDITALIC <- 4;
par(ask=FALSE);

plotYvsFragmentLength <- function(l,y, xlim=c(0,2400), ylim=NULL, xlab="PCR fragment length (bases)", ylab=NULL, main=NULL, ...) {
  # Filter out NAs
  keep <- is.finite(l) & is.finite(y);
  l <- l[keep];
  y <- y[keep];

  smoothScatter(l, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main);
  stext(sprintf("%d SNPs", length(y)), side=3, pos=1, line=-1, cex=0.8);
} # plotYvsFragmentLength()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("modelRmaCnPlm", "modelMbeiCnPlm");
pattern <- ",A[+]B";
if (exists("dataSets", mode="list") && inherits(dataSets[[1]], "CnChipEffectSet")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select arrays to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ces <- dataSets[[1]];
choices <- getNames(ces);
if (!exists("arrays"))
  arrays <- NULL;
arrays <- selectMenu(choices=choices, selected=arrays, title="Select the arrays you which to fit (leave all blank to use all)");
if (length(arrays) == 0) {
  arrays <- seq(choices);
} else {
  arrays <- match(arrays, choices);
} 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");
s <- paste("Arrays: ", seqToHumanReadable(arrays), sep="");
s <- GenericSummary(s);
print(s);


res <- list();
for (kk in seq(along=dataSets)) {
  ces <- dataSets[[kk]];
  verbose && enter(verbose, "Chip-effect set ", getName(ces));

  verbose && cat(verbose, "Chip-effect set:");
  verbose && print(verbose, ces);

  # Hard wired for now
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);

  # Get PCR fragment length data for SNPs
  verbose && enter(verbose, "Get PCR fragment length data for SNPs");
  cdf <- getCdf(ces);
  units <- indexOf(cdf, "^SNP");
  si <- getSnpInformation(cdf);
  fl <- getFragmentLengths(si, units=units);
  # Filter out SNP with missing data
  keep <- is.finite(fl);
  fl <- fl[keep];
  units <- units[keep];

  # Density of all PCR fragment lengths
  d <- density(fl, adjust=1/2);
  # Standardize to maxscale
  d$y <- d$y / max(d$y, na.rm=TRUE);

  # Plot data only for a subset of the SNPs
  subset <- 1:25000;
  units <- units[subset];
  fl <- fl[subset];
  verbose && exit(verbose);

  verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
  ce <- getFile(ces, 1);
  map <- getCellMap(ce, units=units, verbose=less(verbose));
  verbose && str(verbose, map);
  verbose && exit(verbose);

  verbose && enter(verbose, "Retrieving reference signals");
  # Get the reference chip effects
  ceR <- getAverage(ces);
  yR <- getDataFlat(ceR, units=map, fields="theta")[,"theta"];
  verbose && str(verbose, yR);
  lcR <- log2(yR);
  ok <- is.finite(lcR);
  fitR <- lowess(fl[ok],lcR[ok]);
  class(fitR) <- "lowess";
  predictR <- function(x, ...) { predict(fitR, newdata=x, ...) }
  verbose && exit(verbose);

  # Get the output path
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell$", "", chipType);
  path <- filePath("pcrFragmentLength", getFullName(ces), chipType);
  mkdirs(path);

  # Get a fixed 'xlim'
  xlim <- NULL;
  maxFl <- max(fl);
  for (lim in c(2400,1200,600)) {
    if (maxFl <= lim)
      xlim <- c(0,lim);
  }

  # Get a fixed 'lthetaLim'
  lthetaLim <- c(0,16);
  if (regexpr("^Mapping50K", chipType) != -1) {
    lthetaLim <- c(8,18);
  } else if (regexpr("^Mapping250K", chipType) != -1) {
    lthetaLim <- c(5,17);
  }

  # For each array...
  for (array in arrays) {
    # Get chip effects
    ce <- getFile(ces, array);
    name <- getName(ce);
    verbose && enter(verbose, sprintf("Array #%d (%s)", array, name));

    verbose && enter(verbose, "Retrieving signals");
    y <- getDataFlat(ce, units=map, fields="theta")[,"theta"];
    verbose && str(verbose, y);
    verbose && exit(verbose);
    if (length(y) != length(yR)) {
      throw("BUG: 'y' and 'yR' are of different lengths: ", length(y), " != ", length(yR));
      next;
    }

    for (norm in c(FALSE, TRUE)) {
      # Normalize for PCR fragment length effects?
      if (norm) {
        verbose && enter(verbose, "Normalizing for PCR fragment length effects");
        # This requires that 'y' is already available
        lyN <- normalizeFragmentLength(log2(y), fl, targetFcn=predictR);
        yN <- 2^lyN;
        y <- yN;
        verbose && exit(verbose);
      }

      for (doRatio in c(FALSE, TRUE)) {
        tags <- NULL;
        if (norm) {
          tags <- c(tags, "FLN");
        }
        if (doRatio) {
          tags <- c(tags, "CnVsFragLength");
        } else {
          tags <- c(tags, "ThetaVsFragLength");
        }
  
        tags <- paste(tags, collapse=",");
        imgName <- paste(name, tags, sep=",");
        filename <- sprintf("%s.%s", imgName, ext);
        pathname <- filePath(path, filename);
        if (isFile(pathname) && !is.null(dev)) {
          verbose && cat(verbose, "Already plotted...skipping");
          next;
        }
  
        verbose && cat(verbose, "Pathname: ", pathname);

        if (doRatio) {
          lc <- log2(y/yR);
          ylab <- expression(log[2](theta/theta[R]));
          ylim <- c(-1,1)*4;
        } else {
          lc <- log2(y);
          ylab <- expression(log[2](theta));
          ylim <- lthetaLim;
        }

        width <- 1.618*height / 1.08;
        if (!is.null(dev))
          dev(pathname, width=width, height=height);
        opar <- par(font=BOLD, font.axis=BOLD, font.lab=BOLD, font.main=BOLD, font.sub=NORMAL, cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1, mar=c(5,5,2,2)+0.2, lwd=2, pch=1, ps=12, mgp=c(3, 1, 0), tcl=0.3);
  
        plotYvsFragmentLength(fl, lc, xlim=xlim, ylim=ylim, ylab=ylab);
        stext(sprintf("Sample: %s", name), side=3, pos=0, line=0, cex=0.8);
        stextChipType(cdf, line=0);
        if (norm) {
          stext("Normalized for PCR fragment length", side=3, pos=0, line=-1, cex=0.8);
        }

        # Add lowess curves
        ok <- is.finite(lc) & is.finite(fl);
        fit <- lowess(fl[ok],lc[ok]);
        lines(fit, lwd=4, col="red");
  
        if (doRatio) {
          abline(h=0, lwd=4, col="purple");
        } else {
          lines(fitR, lwd=4, col="purple");
        }
  
        # Add the density plot of fragment lengths
        ylim <- par("usr")[3:4];
        dd <- d; dd$y <- ylim[1] + dd$y * (diff(ylim) / 10);
        lines(dd, col="darkgray", lwd=2);
    
        par(opar);
        if (!is.null(dev))
          dev.off();
      } # for (doRatio in ...)
    } # for (norm in ...)
   
    verbose && exit(verbose);
  } # for (array ...)
}
verbose && exit(verbose);


############################################################################
# HISTORY:
# 2006-11-28
# o Recreated.
############################################################################
