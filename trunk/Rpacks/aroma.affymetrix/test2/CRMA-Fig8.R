# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (interactive())
  savehistory();
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");
log <- Verbose(threshold=-5, timestamp=TRUE);


imgFormat <- "screen";
imgFormat <- "eps";
force <- FALSE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getBlockAverageMap <- function(n, h=1, s=0, ...) {
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(1,1000));


  # Is h an integer?
  if (h == as.integer(h)) {
    idxs <- matrix(seq_len(n), nrow=h);
    if (n %% h != 0)
      idxs <- idxs[,-ncol(idxs)];
  } else {
    K <- ceiling(n/h);
    idxs <- matrix(TRUE, nrow=ceiling(h), ncol=K);
    steps <- (h %% 1) * ceiling(K);
    incl <- seq(from=1, to=K, length=steps);
    incl <- round(incl);
    idxs[ceiling(h), -incl] <- FALSE;
    # Shift
    if (s > 0) {
      lastRow <- idxs[ceiling(h),];
      tail <- seq(from=length(lastRow)-s+1, to=length(lastRow));
      lastRow <- c(lastRow[tail], lastRow[-tail]);
      idxs[ceiling(h),] <- lastRow;
    }
    idxs[idxs] <- seq_len(n);
    idxs[idxs == 0] <- NA;
  }

  # Skip last column in case looping
  if (n %% h != 0)
    idxs <- idxs[,-ncol(idxs)];

  # The effective 'h'
  hApprox <- sum(!is.na(idxs))/ncol(idxs);
  attr(idxs, "hApprox") <- hApprox;

  idxs;
} # getBlockAverageMap()


blockAvg <- function(X, idxs, FUN=rowMedians.matrix, ...) {
  na.rm <- (any(is.na(X)) || any(is.na(idxs)));
  dimnames <- dimnames(X);
  dimnames(X) <- NULL;
  X <- t(X);
  Z <- apply(idxs, MARGIN=2, FUN=function(jj) {
    jj <- jj[is.finite(jj)];
    Zjj <- X[,jj,drop=FALSE];
    Zjj <- FUN(Zjj, na.rm=na.rm);
    Zjj;        
  });
  Z <- t(Z);
  colnames(Z) <- dimnames[[2]];
  Z;
} # blockAvg()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetName <- "Affymetrix_2006-HapMap270.CEU.founders";
dataPath <- filePath("data", dataSetName);

if (!exists("truth")) {
  filename <- "etc.RData";
  pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                     mustExist=TRUE);
  attachLocally(loadObject(pathname));

  filename <- "truth.RData";
  pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                     mustExist=TRUE);
  attachLocally(loadObject(pathname));

#  cdf <- AffymetrixCdfFile$fromChipType("Mapping250K_Nsp");
#  gs <- getUnitSizes(cdf, units=units);
}

males <- which(truth[1,] == 0);
females <- which(truth[1,] == 1);


# Subset of arrays to analyze
excl <- c("NA12145");
cc <- cc0 <- which(!colnames(truth) %in% excl);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

R <- 5;

case <- "optimal";

# Classification
group <- c("male"=0, "female"=1);
toCall <- 0;
pName <- names(group[group == toCall]);
nName <- names(group[group != toCall]);
tpName <- pName;
ylab <- sprintf("True-positive rate\n(correctly calling %ss %ss)", pName, pName);
xlab <- sprintf("False-positive rate\n(incorrectly calling %ss %ss)", nName, pName);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
useColors <- FALSE;
# Setup a palette
if (useColors) {
  niceCols <- RColorBrewer::brewer.pal(12, "Paired");
  niceCols <- niceCols[c(6,2,8,4)];
  useColorsStr <- "-col";
} else {
  niceCols <- c("#000000", "#444444", "#888888", "#aaaaaa");
  useColorsStr <- "";
}

fig <- 1;


# Arrays with strong spatial artifacts (in residuals)
badArrays <- c("NA06985", "NA07055", "NA11829", "NA11830", "NA11995", "NA12004", "NA12005", "NA12006", "NA12057", "NA12234", "NA12156", "NA12236", "NA12239", "NA12264", "NA12760", "NA12762", "NA12763", "NA12892");

if (!exists("sets", mode="list")) {
  sets <- list(
   "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[1], name="CRMA"),
   "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
   "CNAG,FLN+GC,r60,xport" = list(col=niceCols[3], name="CNAG*"),
   "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT")
  );
}



references <- seq_len(ncol(truth));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=sets)) {
  name <- names(sets)[kk];
  verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", 
                                                   kk, name, length(sets)));

  # Get set
  set <- sets[[name]];
  if (is.null(set$M)) {
    # Missing 'theta'?
    if (is.null(set$theta)) {
      # Load thetas for chromosome X
      dataName <- gsub(";.*$", "", name);
      filename <- sprintf("%s.RData", dataName);
      pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                             mustExist=TRUE);
      print(filename);
      data <- loadObject(pathname);
      if ("theta" %in% names(data))
        set$theta <- data$theta;
      if ("M" %in% names(data))
        set$M <- data$M;
      rm(data);
    }

    # Has 'theta'?
    if (!is.null(set$theta)) {
      # Calculate average theta, and (A,M).
      thetaRefs <- set$theta[,references];
      set$thetaR <- apply(thetaRefs, MARGIN=1, FUN=median, na.rm=TRUE);
      str(thetaRefs);

      set$M <- log2(set$theta/set$thetaR); 
      set$A <- log2(set$theta*set$thetaR)/2;
    } else if (!is.null(set$M)) {
      if (name == "CNAT,QN,PLIER,FLN+GC,xport") {
        # Make sure order of samples are the same as in 'truth'
        df <- read.table("CNAT/foo.xls", header=TRUE, sep="\t");
        expId <- sprintf("%s.cn", df[,1]);
        sampleNames <- gsub("[.]cel$", "", df[,2]);
        cc <- match(colnames(set$M), expId);
        stopifnot(all(is.finite(cc)));
        sampleNames <- sampleNames[cc];
        cc <- match(sampleNames, colnames(truth));
        stopifnot(all(is.finite(cc)));
        set$M <- set$M[,cc];
        rm(df, expId, sampleNames, cc);
      }
    }

    # Assert that we at least got an 'M' field
    stopifnot("M" %in% names(set));

    sets[[name]] <- set;
  }

  verbose && exit(verbose);
} # for (kk in ...)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Average in non-overlapping blocks 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
keep <- (x > 2.8e6);
cc <- cc0;

fig <- 1;

if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig);

  swapXY <- FALSE;
  #swapXY <- TRUE;
  hs <- seq(from=1, to=7.5, by=1/4);
#  hs <- seq(from=1, to=7.5, by=1);
#  hs <- c(hs, 1.85, 1.94, 1.95, 1.96, 3.706, 5.555, 7.43);
  hs <- sort(hs);
  xMean <- mean(diff(sort(x[keep])));

  xlim <- c(xMean, max(xMean*hs))/1e3;
  xlim <- c(xMean, 200e3)/1e3;
#  ylim <- c(0.6,1);
  ylim <- c(0.3,1);
#  ylim <- c(0,1);
  fpRate <- 1/sum(keep);
  fpRate <- 1/2*1/29;
  fpRate <- 1/100;
  
#  imgFormat <- "screen";
  
  spar <- 0.48;
#  spar <- 0.55;
#  spar <- 0.65;
#  spar <- 0.2;
  spar <- 0.35;
  xlab <- "Average distance between window centers (kb)";
#  ylab <- "TP rate";
  xSeq <- seq(0, xlim[2], by=10);
  xSeq <- xSeq[xSeq > xMean/1e3];
  if (length(xSeq) > 100) stop("!!!");
  
  if (imgFormat == "eps") {
    imgName <- sprintf("TPvResolution,%s,fpRate%.4f", case, fpRate);
    imgName <- gsub("[.]", "_", imgName);
    filename <- sprintf("%s%s.%s", imgName, useColorsStr, imgFormat);
    eps(filename, width=6, height=1.33*0.618*6);
    par(mar=c(2.8,3.8,3,0.8)+0.1);
    par(cex=1.2, mgp=c(1.7,0.5,0));
    par(cex.axis=1, cex.lab=1, xaxs="i", yaxs="i", xpd=FALSE);
  } else if (imgFormat == "screen") {
    par(mar=c(4,4.5,3.5,1)+0.1);
    par(cex=1, mgp=c(2.5,0.7,0));
    par(cex.axis=0.8, xaxs="i", yaxs="i", xpd=FALSE);
  }
  
  
  if (swapXY) {
    plot(NA, xlim=ylim, ylim=xlim, ylab=xlab, xlab=ylab, axes=FALSE);
    abline(v=seq(ylim[1], ylim[2], by=0.1), lwd=1, lty=3);
    abline(h=xSeq, lwd=1, lty=3);
  } else {
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
    abline(h=seq(ylim[1], ylim[2], by=0.1), lwd=1, lty=3);
    abline(v=xSeq, lwd=1, lty=3);
#    abline(v=27*1:7, lwd=1, lty=3);
  }
  axis(side=1);  
#   axis(side=1, at=27*1:7, labels=27*1:7);
   axis(side=2);
  stext(side=4, pos=0, sprintf("FP rate: %.2f%%", 100*fpRate), cex=0.7);
  if (swapXY) {
    abline(h=xMean/1e3, lwd=2, lty=1);
    axis(side=2, at=xMean/1e3, sprintf("%.1f", xMean/1e3));
    abline(v=ylim);
    axis(side=4, at=xMean/1e3*hss, hss);
  } else {
#    abline(v=xMean/1e3, lwd=2, lty=1);
#    abline(h=ylim);
    axis(side=1, at=xMean/1e3, sprintf("%.1f", xMean/1e3));
    hss <- seq(from=1, to=ceiling(max(hs)));
    axis(side=3, at=xMean/1e3*hss, hss);
    stext(side=3, pos=0.5, "Average number of loci per window", line=1.5, cex=1.2);
  }
  
  par(xpd=FALSE);
  ltys <- rep(c(1,6), length.out=length(sets));
  for (kk in seq(along=sets)) {
    name <- names(sets)[kk];
    if (regexpr("R=[0-9]", name) != -1)
      next;
  
    verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", 
                                                     kk, name, length(sets)));
    set <- sets[[name]];
    M <- set$M;
    hTpName <- sprintf("hTp_%.4f", fpRate);
    hTp <- set[[hTpName]];
    if (force)
      hTp <- NULL;
    hTp <- scanTpAtFp(truth[keep,cc], M[keep,cc], x[keep], hs=hs, hTp=hTp, fpRate=fpRate, call=toCall, shifts=seq(from=0, to=(max(hss)-1), by=3), verbose=verbose);
    set[[hTpName]] <- hTp;
  
    hTp[,c("xMean", "xMedian")] <- hTp[,c("xMean", "xMedian")]/1e3;
  
    xField <- "xMean";
  #  xField <- "xMedian";
    fit <- smooth.spline(hTp[,xField], hTp[,"tpRate"], all.knots=FALSE, spar=spar);
  #  fit <- loess(hTp[,"tpRate"] ~ hTp[,xField], span=0.3);
  #  fit$y <- fit$fitted;
    fit$y[fit$y > 1] <- 1;
    fit$y[fit$y < 0] <- 0;
  
#    points(hTp[,xField], hTp[,"tpRate"], pch=19, col=set$col, cex=1);
  
    if (spar < 0.37) {
      par(xpd=TRUE);
#      points(hTp[1,xField], hTp[1,"tpRate"], pch=19, col=set$col, cex=1.3);
      par(xpd=FALSE);
    }
    if (swapXY) {
      lines(fit$y, fit$x, lwd=3, col=set$col, lty=ltys[kk]);
    } else {
      lines(fit$x, fit$y, lwd=3, col=set$col, lty=ltys[kk]);
    }
  
    sets[[name]] <- set;
  
    verbose && exit(verbose);
  } # for (kk in ...)
  
  cols <- sapply(sets, FUN=function(set) set$col);
#  ltys <- sapply(sets, FUN=function(set) set$lty);
#  ltys <- 1;
  names <- sapply(sets, FUN=function(set) set$name);
  labels <- names;
  legend("bottomright", lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");
  
  par(xpd=FALSE);
  if (swapXY) {
    abline(v=ylim);
  } else {
    abline(h=ylim);
  }
  
  if (imgFormat != "screen") {
    dev.off();
  }
}
