library(R.oo);

myPlotDensity <- function(cs, ...) {
  if (imgFormat == "eps") {
    imgName <- sprintf("%s,plotDensity", getFullName(cs));
    filename <- sprintf("%s.eps", imgName);
    filename <- gsub("[.]", "_", filename);
    eps(filename, width=5, height=0.618*5);
    par(mar=c(4,5,1,1)+0.1);
    par(cex=0.7, mgp=c(2.7,0.7,0));
    on.exit(dev.off());
  } else if (imgFormat == "screen") {
    Device$set(2, height="61.8%");
    Device$setStyle("PowerPoint");
    par(mar=c(5,5,2,1)+0.1);
    par(cex=1, mgp=c(2.5,0.2,0));
  }
  par(cex.axis=1, cex.lab=1, font.lab=2);
  plotDensity(cs, subset=1/5, types="pm", ylim=c(0,0.4), annotate=FALSE, verbose=verbose);
} # myPlotDensity()


setMethodS3("plotAllelicCrosstalk", "AffymetrixCelFile", function(this, pch=NA, lwd=4, col="red", xlim=c(-2^9,2^12), ylim=xlim, xaxs="r", yaxs=xaxs, transformation=function(x) x^0.18, imgFormat=c("screen", "eps"), ...) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  myPlotAllelicCrosstalk <- function(y, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, basepair=c("A", "B"), tags=NULL, ...) {
    if (is.null(xlab)) {
      xlab <- substitute(y[A], list=list(A=basepair[1]));
      xlab <- substitute(PM[A], list=list(A=basepair[1]));
    }
    if (is.null(ylab)) {
      ylab <- substitute(y[B], list=list(B=basepair[2]));
      ylab <- substitute(PM[B], list=list(B=basepair[2]));
    }
  
    y[y < xlim[1]-0.1*diff(xlim) | y > xlim[2] + 0.1*diff(xlim)] <- NA;


    if (imgFormat == "eps") {
      tags <- paste(c("", tags), collapse=",");
      imgName <- sprintf("%s,ACC,%s%s%s", getFullName(this), basepair[1], basepair[2], tags);
      imgName <- gsub("[.]", "_", imgName);
      imgName <- gsub(",,", ",", imgName);
      filename <- sprintf("%s.eps", imgName);
      eps(filename, width=5, height=5);
      par(mar=c(4,5,1,1)+0.1);
      par(cex=0.7, mgp=c(2.7,0.7,0));
    } else if (imgFormat == "screen") {
      Device$set(2, height="61.8%");
      Device$setStyle("PowerPoint");
      par(mar=c(4,5,1,1)+0.1);
      par(cex=1);
      par(cex.lab=1.8, mgp=c(2.5,0.2,0));
    }

    smoothScatter(y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...);
  } # myPlotAllelicCrosstalk()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'imgFormat':
  imgFormat <- match.arg(imgFormat);  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add plot hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  myLinesAllelicCrosstalk <- function(a, B, lwd=4, col="red", max=1e5, ...) {
    print(a);
    print(B);
    bA <- B[1,2]/B[1,1]; 
    bB <- B[2,1]/B[2,2]; 
    lines(x=a[1]+c(0,1)*max, y=a[2]+c(0,1)*max*bA, col=col, lwd=lwd, ...);
    lines(x=a[1]+c(0,1)*max*bB, y=a[2]+c(0,1)*max, col=col, lwd=lwd, ...);
  } # myLinesAllelicCrosstalk()


  hookName <- "calibrateAllelicCrosstalk.AffymetrixCelFile.onData";
  setHook(hookName, function(this, y, basepair, ...) { 
    myPlotAllelicCrosstalk(y=y, pch=pch, xlim=xlim, ylim=ylim, basepair=basepair, transformation=transformation, xaxs=xaxs, yaxs=yaxs); 
  }, action="replace");

  hookName <- "calibrateAllelicCrosstalk.AffymetrixCelFile.onFit";
  setHook(hookName, function(this, fit, ...) { 
    myLinesAllelicCrosstalk(a=fit$origin, B=fit$W, lwd=lwd, col=col);
    if (imgFormat != "screen")
      dev.off();
  }, action="replace");

  hookName <- "calibrateAllelicCrosstalk.AffymetrixCelFile.onUpdated";
  setHook(hookName, function(this, y, basepair, fit, yC, ...) { 
    myPlotAllelicCrosstalk(y=yC, pch=pch, xlim=xlim, ylim=ylim, basepair=basepair, transformation=transformation, xaxs=xaxs, yaxs=yaxs, tags="ACC"); 
    myLinesAllelicCrosstalk(a=c(0,0), B=matrix(c(1,0,0,1),ncol=2), lwd=lwd, col=col);
    if (imgFormat != "screen")
      dev.off();
  }, action="replace");


  path <- tempdir();
  dfC <- calibrateAllelicCrosstalk(this, path=path, verbose=verbose);
  pathname <- getPathname(dfC);
  file.remove(pathname);
}, protected=TRUE);

stop()

library(aroma.affymetrix);
source("../test/init.R");
verbose <- Arguments$getVerbose(-5);
timestampOn(verbose);



imgFormat <- "screen";
imgFormat <- "eps";

fig <- 1;

name <- "Affymetrix_2006-HapMap270.CEU.founders";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

for (chipType in chipTypes[1]) {
  print(chipType);

  if (!exists("cs")) {
    path <- filePath("rawData", name, chipType, expandLinks="any");
    cs0 <- AffymetrixCelSet$fromFiles(path, verbose=verbose);
  
    # Get all parents in the CEU population
    population <- sapply(cs0, FUN=getAttribute, "population");
    fatherID <- as.integer(sapply(cs0, FUN=getAttribute, "fatherID"));
    idxs <- which((population == "CEU") & is.na(fatherID));

print(idxs);  
    cs <- cs0;
#    cs <- extract(cs0, idxs);
#    setAlias(cs, sprintf("%s.CEU.founders", getName(cs0)));
    rm(cs0);
    n23 <- as.integer(sapply(cs, FUN=getAttribute, "n23"));
  }
  print(cs);

stop()

# Plot density  
#    myPlotDensity(cs);


  cdf <- getCdf(cs);
  gi <- getGenomeInformation(cdf);
  units <- getUnitsOnChromosome(gi, 23);

  subsetToFit <- seq_len(nbrOfUnits(cdf));
  subsetToFit <- setdiff(subsetToFit, units);

#  csC <- process(acc, verbose=verbose);
#  bg <- RmaBackgroundCorrection(cs);

  acc <- AllelicCrosstalkCalibration(cs);
  print(acc);
  csC <- process(acc, verbose=verbose);
  print(csC);

#  qn <- DChipQuantileNormalization(csC, typesToUpdate="pm", 
#                                                 tags=c("*", "pm", "-X"));
#  excludeChrXFromFit(qn);
#  print(qn);
#  csN <- process(qn, verbose=verbose);
#  print(csN);

  gc();

#  plm <- MbeiCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("-100","*"));
  plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, tags=c("+300","*"));
  plm$treatNAsAs <- "NA";
  plm$shift <- +300;
  print(plm);
  ces <- getChipEffectSet(plm);
  fit(plm, ram=1.5, verbose=verbose);

#  fit(plm, units=units, ram=2.5, verbose=verbose);

  fln <- FragmentLengthNormalization(ces, subsetToFit=subsetToFit, tags=c("*", "-X"));
  print(fln);
  cesN <- process(fln, verbose=verbose);
  print(cesN);

#  gcn <- GcContentNormalization(cesN, subsetToFit=subsetToFit, tags=c("*", "-X"));
#  cesN2 <- process(gcn, verbose=verbose);

  ces <- cesN;

stop();

  qn2 <- DChipQuantileNormalization(cesN, tags=c("*", "-X"));
  excludeChrXFromFit(qn2);
  print(qn2);

  cesQN <- process(qn2, verbose=verbose);
#  print(csN);

stop()
ces <- cesQN;
  print(ces);
  if (TRUE) {
    rm(theta); 
    theta <- extractMatrix(ces, units=units, verbose=verbose);
    saveObject(list(theta=theta), sprintf("%s.RData", paste(getTags(ces), collapse=",")), compress=FALSE);
    rm(theta); 
  }

  if (!plm$mergeStrands) {
    merger <- ChipEffectGroupMerge(cesN);
    setMergeFunction(merger, function(x) { 
      x[1,] <- 2^colMeans(log2(x), na.rm=TRUE); 
      x[-1,] <- 0;
      x;
    })
    cesM <- process(merger, verbose=verbose);
  }

  gc();
}
