# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# First some utilities functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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


plotPairDensities <- function(rg, xlim=c(-2,18), ylim=c(0,0.8), ...) {
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
  if (!is.null(names)) {
    usr <- par("usr")
    legend(x=usr[2],y=usr[4], legend=names, fill=colors, xjust=1, yjust=1, cex=0.7);
  }
  par(opar);
}




# Draw the (R,G) grid of the (fitted) affine model
drawRGGrid <- function(maxSignal=16, by=1, drawCurve=TRUE, highlightLog1=FALSE, afit=NULL, aR=NULL, aG=NULL, b=NULL) {
  x <- seq(-2*maxSignal,maxSignal, by=by)
  # The grid ticks on the non-logarithmic scale and shifted -1.
  X <- 2^x
  # The input (R,G) grid.
  R <- matrix(X, nrow=length(X), ncol=length(X))
  G <- t(R)

  if (!is.null(afit)) {
    if (is.null(aR))
      aR <- afit$a[1];
    if (is.null(aG))
      aG <- afit$a[2];
    if (is.null(b))
      b <- max(afit$b[-1]);
  }

  if (!is.null(aR) & !is.null(aG) & !is.null(b)) {
    R <- aR +   R
    G <- aG + b*G
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
      lines(a[z,],m[z,], col=col[kk], lwd=lwd[kk]);
      lines(a[,z],m[,z], col=col[kk], lwd=lwd[kk]);
    }
  }

  if (drawCurve) {
    lines(diag(a), diag(m), col="blue", lty=4, lwd=2);
  }
} # drawRGGrid()

