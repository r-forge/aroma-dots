plotPairwiseScans <- function(rg, channel, ...) {
  pair <- NULL;
  for (ii in 1:(nbrOfSlides(rg)-1)) {
    for (jj in (ii+1):nbrOfSlides(rg)) {
      pair <- cbind(pair, c(ii,jj));
    }
  }

  opar <- par(mar=c(5,5,3,1));
  xlab <- substitute(y[channel]^{(v)}, list=list(channel=channel));
  ylab <- substitute(y[channel]^{(w)}, list=list(channel=channel));
  plot(NA, xlab=xlab, ylab=ylab, ...);
  box(lwd=2, col=switch(channel, "R"="red", "G"="green"));
  colors <- 1:ncol(pair)+3;
  for (kk in 1:ncol(pair)) {
    s <- pair[1,kk];
    t <- pair[2,kk];
    yv <- rg[[channel]][,s];
    yw <- rg[[channel]][,t];
    ok <- (yv > 0 & yw > 0 & is.finite(yv) & is.finite(yw));
    mc <- log((yv/yw)[ok], base=2);
    ac <- log((yv*yw)[ok], base=2)/2;
    points(ac,mc, pch=".", col=colors[kk]);
  }

  # Add a legend
  names <- getSlideName(rg);
  pairNames <- apply(pair, MARGIN=2, FUN=function(x) paste(names[x], collapse=","));
  pairNames <- paste("(", pairNames, ")", sep="");
  usr <- par("usr")
  legend(x=usr[2],y=usr[3], legend=pairNames, fill=colors, xjust=1, yjust=0, cex=0.7);
  par(opar);
}


plotPairDensities <- function(rg, xlim=c(2,18), ylim=c(0,0.8), ...) {
  colors <- seq(rg)+3;
  opar <- par(mar=c(5,5,3,1));
  xlab <- expression(log[2](y[c]));
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab="density", ...);
  for (ch in c("R", "G")) {
    col <- switch(ch, "R"="red", "G"="green");
    for (slide in seq(rg)) {
      x <- rg[[ch]][,slide];
      line <- density(log(x[is.finite(x) & x > 0], base=2));
      lines(line, col=col, lwd=4);
      lines(line, col=colors[slide], lwd=1);
    }
  }

  # Add a legend
  names <- getSlideName(rg);
  pairNames <- paste(names, "V", sep="");
  usr <- par("usr")
  legend(x=usr[2],y=usr[4], legend=pairNames, fill=colors, xjust=1, yjust=1, cex=0.7);
  par(opar);
}



# Draw the (R,G) grid of the (fitted) affine model
drawRGGrid <- function(maxSignal=16, by=1, drawCurve=TRUE, highlightLog1=TRUE, afit=NULL) {
  x <- seq(-maxSignal,maxSignal, by=by)
  # The grid ticks on the non-logarithmic scale and shifted -1.
  X <- 2^x
  # The input (R,G) grid.
  R <- matrix(X, nrow=length(X), ncol=length(X))
  G <- t(R)

  if (!is.null(afit)) {
    a <- getADiagonal(afit)[1];
    R <- a + R
    G <- a + G
  }

  r <- log(R, base=2)
  g <- log(G, base=2)
  m <- r-g
  a <- (r+g)/2
  drawGrid(a,m, col="gray");
  if (highlightLog1) {
    z <- which(x == 0);
    col <- c("lightblue", "gray");
    lwd <- c(2,1);
    for (kk in 1:2) {
      lines(a[z,],m[,z], col=col[kk], lwd=lwd[kk]);
      lines(a[,z],m[z,], col=col[kk], lwd=lwd[kk]);
    }
  }

  if (drawCurve)
    lines(diag(a), diag(m), col="red")
} # drawRGGrid()

arrays <- c("12633506", "12633508", "12638525", "12638526", "12638528", "12638529");
scanners <- c("Agilent", "AxonGenePix4000A");
imethods <- c("gpr", "spot");
for (array in arrays) {
  for (scanner in scanners) {
    for (imethod in imethods) {
rm(rg);

# One array was scanned four times at four different PMT settings.
if (!exists("rg")) {
  if (array %in% c("12633506", "12633508"))
    date <- "20031007"
  else
    date <- "20031015";
  path <- file.path(getOption("DATAROOT"), "Oncology", date, scanner, "non-calibrated");
  filename <- paste(array, imethod, "RawData", "dat", sep=".");
  print(paste(path, filename, sep="/"));
  rg <- RGData$read(filename, path=path);
#  subset <- seq(1,nbrOfSpots(rg),length=5000);  keepSpots(rg, subset);
  setSlideName(rg, c("500V","600V","700V","800V"));
  #rg <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
  #setLayout(rg, Layout(4,4,17,17));       # Not really necessary!
  rm(fits);
}

# Reorder the slides in scan order
keepSlides(rg, c("800V","500V","600V","700V"));

# Remove suspicious slides
#removeSlides(rg, "700V");

# Number of scans
K <- nbrOfSlides(rg);

# Number of bootstrap samples
R <- 100;

# Bootstrap quantiles to be highlighted
probs <- c(0.9,0.95,0.99);
probsStr <- sapply(probs, FUN=sprintf, fmt="%.2f");


timestamp <- format(Sys.time(), format="%Y%m%d-%H%M%S");
resfile <- sprintf("bootstrap-%s.dat", timestamp);

Device$setPrintPath(file.path("figures", scanner));
fig <- 1;
channels <- c("R", "G");
if (!exists("fits")) {
  fits <- fitMultiscanAffine(rg, bootstrap=TRUE, R=R, verbose=TRUE);
}
for (channel in channels) {
  fig <- fig + 1;
  Device$set(fig);
  col <- c(R="red", G="green")[channel];

  # Get initial estimate of \gamma_0 = (a_0 - \alpha_0)
  fit <- fits[[channel]];
  print(fit);
  df <- fit$t;
  df <- cbind(channel, df);
  df <- cbind(imethod, df);
  df <- cbind(scanner, df);
  df <- cbind(array, df);
  df <- cbind(df, nbrOfSpots=nbrOfSpots(rg));
  append <- file.exists(resfile);
  write.table(file=resfile, df, sep="\t", quote=FALSE, row.names=FALSE, col.names=!append, append=append);

  a <- fit$t[,(1:K)];
  alpha <- fit$t[,(1:K)+2*K];
  d <- a-alpha;
  colnames(a) <- paste("a", getSlideNames(rg));
  colnames(alpha) <- paste("alpha", getSlideNames(rg));
  colnames(d) <- paste("d", getSlideNames(rg));

  pair <- NULL;
  for (ii in 1:(nbrOfSlides(rg)-1)) {
    for (jj in (ii+1):nbrOfSlides(rg)) {
      pair <- cbind(pair, c(ii,jj));
    }
  }

  xlim <- ylim <- c(10,30);
  nbrOfPairs <- ncol(pair);

  subplots(nbrOfPairs+2);
  opar <- par(mar=c(5,5,3,1));

  for (kk in 1:nbrOfPairs) {
    s <- pair[1,kk];
    t <- pair[2,kk];
    sname <- getSlideNames(rg)[s];
    tname <- getSlideNames(rg)[t];
    yv <- rg[[channel]][,s];
    yw <- rg[[channel]][,t];
    ok <- (is.finite(yv) & is.finite(yw));

    # Plot the data
    xlab <- substitute(y[channel]^{(v)}, list=list(channel=channel, v=sname));
    ylab <- substitute(y[channel]^{(w)}, list=list(channel=channel, w=tname));
    plot(yv[ok],yw[ok], pch=".", col=col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
    box(col=col);
    abline(a=0,b=1, col="lightgray")

    # Plot the a estimates
    points(a[,s],a[,t], pch=".", col="black");
#    points(mean(a[,s]),mean(a[,t]), pch="o", col="black", lwd=2, cex=1.5);

    # Plot the alpha estimates
    points(alpha[,s],alpha[,t], pch=".", col="blue");
#    points(mean(alpha[,s]),mean(alpha[,t]), pch="+", col="blue", lwd=2, cex=1.5);
  }

  # Boxplot of the estimated a vector
  boxplot(as.data.frame(a), las=2, ylim=ylim);
  box(col=col);

  boxplot(as.data.frame(d), las=2, ylim=c(-1,1)*3.5);
  box(col=col);

  # Parallel plot of the relative a vector
  o <- order(d[,1]);
  bcol <- getColors(TopologyColor(seq(0,1,length=nrow(a))));
  plotParallel(d[o,], col=bcol);
  box(col=col);

  imgname <- sprintf("%s.%s-Y%svsY%s-a-alpha-zoom.eps", array, imethod, channel, channel);
  Device$print(imgname, width="150%", height="160%");
} # for (channel ...)

    } # for (imethod in imethods) 
  } # for (scanner in scanners)
} # for (array in arrays)
  

