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


# One array was scanned four times at four different PMT settings.
rg <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
setLayout(rg, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg, c("500V","600V","700V","800V"));
# Reorder the slides in scan order
keepSlides(rg, c("800V","500V","600V","700V"));

# Remove suspicious slides
# removeSlides(rg, "500V");

# Number of scans
K <- nbrOfSlides(rg);

# Number of bootstrap samples
R <- 200;

# Bootstrap quantiles to be highlighted
probs <- c(0.9,0.95,0.99);
probsStr <- sapply(probs, FUN=sprintf, fmt="%.2f");

channels <- c("R", "G");
subplots(length(channels), nrow=2);
for (channel in channels) {
  col <- c(R="red", G="green")[channel];

  # Get initial estimate of \gamma_0 = (a_0 - \alpha_0)
  fit0 <- fitMultiscanAffine(rg, channel=channel)[[1]];
  print(fit0);
  # gamma <- (a - alpha)
  gamma0 <- getDistance(fit0);
  print(gamma0);
  T0 <- sqrt(sum(gamma0^2))
 
  # Generate sample under H0 using bootstrap
  # Shift by -gamma0 == -(a-alpha) == (alpha-a). So if data is under H1 we sample
  # from shifted y_i = a + z_i + (alpha-a) == alpha + z_i, i.e. under H0.
  fitH0 <- fitMultiscanAffine(rg, channel=channel, aShift=-gamma0, bootstrap=TRUE, R=R)[[1]];
  print(fitH0);
  aH0 <- fitH0$t[,(1:K)];
  alphaH0 <- fitH0$t[,(1:K)+2*K];
  gammaH0 <- aH0-alphaH0; colnames(gammaH0) <- NULL;
  TH0 <- sqrt(rowSums(gammaH0^2))

  # Plot the bootstrap density
  d <- density(TH0);
  xlim <- range(c(d$x,T0), na.rm=TRUE);
  plot(d, lwd=2, col=col, xlim=xlim, xlab=expression(T==abs(a-alpha)), ylab="density", main=paste("Statistica for channel", channel));

  # Plot the quantiles
  ymid <- mean(par("usr")[3:4]);
  Tq <- quantile(TH0, probs=probs);
  for (kk in 1:length(probs)) {
    abline(v=Tq[kk], lwd=1, lty=3, col=col);
    label <- substitute(T[p] == Tq, list=list(p=probsStr[kk], Tq=sprintf("%.2f", Tq[kk])));
    text(x=Tq[kk], y=ymid, label=label, cex=0.8, srt=90, adj=c(0,-0.5));
  }

  # Plot the estimate
  abline(v=T0, lwd=2, col=col);
  label <- substitute(T[0] == T0, list=list(T0=sprintf("%.2f", T0)));
  text(x=T0, y=ymid, label=label, cex=0.8, srt=90, adj=c(0,-0.5));

} # for (channel ...)
  

