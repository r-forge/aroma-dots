%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{calibrateMultiscan.RGData}
\alias{calibrateMultiscan.RGData}
\alias{RGData.calibrateMultiscan}
\alias{calibrateMultiscan.RGData}
\alias{calibrateMultiscan,RGData-method}


\title{Calibrates multiple re-scanned images based on an affine model}

\description{
  Calibrates multiple re-scanned images based on an affine model.

  Each channel is calibrated seperately.
}

\usage{\method{calibrateMultiscan}{RGData}(this, slides=NULL, channels=NULL, groupBy=NULL, ...)}

\arguments{
 \item{slides}{Slides to be used in the fit \emph{and} that are to be
  calibrated. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{channels}{\code{\link[base]{character}} string specifying which channels to be
   calibrated. If \code{\link[base]{NULL}}, all channels are calibrated.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link[aroma]{LayoutGroups}}
    specifying the groups of spots that are to calibrated individually.
    If \code{\link[base]{NULL}}, all spots are considered.}
 \item{...}{Additional arguments accepted by
    \code{\link[aroma.light]{calibrateMultiscan.matrix}}.
    Its help page is very USEFUL.}
}

\value{
  Returns a \code{\link[base]{list}} containing one element for each calibrated channel.
  Each channel element contains parameter estimates either directly
  as global estimates or as a \code{\link[base]{list}} consisting estimates for each
  group (as defined \code{groupBy}).
}

\section{Weighted calibration}{
 Each data point can be assigned a weight in [0,1] specifying how much
 it should \emph{affect the fitting of the calibration function}.
 Note that here a \emph{data point} is here considered to be the \code{\link[base]{vector}}
 of values from all scans ("slides").

 Regardless of weights, \emph{all} data points are \emph{calibrated} based
 on the fitted calibrated function.

 Data-point weights are obtained from probe weights, if given.
 Weights can be set using \code{\link[aroma:setProbeWeights.RGData]{*setProbeWeights}()}.
 If weights are specified, they will be used.
 Currently it is not possible to set different in different channels.
}

\details{
 Fitting is done by iterated re-weighted principal component analysis
 (IWPCA).
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] H. Bengtsson, J. Vallon-Christersson and G. Jönsson, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr
}

\examples{
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


 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Main example code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# One array was scanned four times at four different PMT settings.
rg0 <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
setLayout(rg0, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg0, c("500V","600V","700V","800V"));
# Reorder the slides (to make the M vs A plot more clear!)
keepSlides(rg0, c("800V","500V","600V","700V"));

for (ch in c("R", "G")) {
  signalWeights <- SignalWeights$fromSaturatedSignals(rg0[[ch]]);
  setSignalWeights(rg0, channel=ch, weights=signalWeights);
  rm(signalWeights);
}

# Calibrate using x <- (y-a)/b. For illustration purposed, we want to
# keep all scans so we do not calculate the average scan here.
rgC <- clone(rg0);
fit <- calibrateMultiscan(rgC, average=NULL, project=FALSE);
str(fit);

Device$set(2, height="108\%");
subplots(9);

Alim <- c(0,20);
Mlim <- c(-1,1)*abs(diff(range(Alim)));

# Same scale on calibrated and non-calibrated data
xlim <- c(4,18); ylim <- c(-1,1)*3.5;

for (ch in c("R", "G")) {
  plotPairwiseScans(rg0, channel=ch, xlim=Alim, ylim=Mlim, main="non-calibrated");
  drawRGGrid(afit=fit[[ch]]);
}

plotPairDensities(rg0, main="non-calibrated");

# It is clear from the calibrated plots that the 500V scan is special.
# By visual inspection of the density plots the following manual
# correction makes their densities "nicer". Uncomment to test.
# rg$R[,1] <- rg$R[,1] - 94; rg$G[,1] <- rg$G[,1] - 54;

for (ch in c("R", "G")) {
  plotPairwiseScans(rgC, channel=ch, xlim=Alim, ylim=Mlim, main="calibrated");
  drawRGGrid(maxSignal=Alim[2]);
}

plotPairDensities(rgC, main="calibrated");


# Now when we have shown the effect of multiscan calibration, we
# calibrate and calculate the median scan as done in practice.
rgCAvg <- clone(rg0);
calibrateMultiscan(rgCAvg);
setSlideNames(rgCAvg, "calib");


# Put both raw and average calibrated data into the same MAData object
rgN <- clone(rgCAvg);
fit <- normalizeAffine(rgN);

# Try also this to see the effect of lowess normalization
maN <- as.MAData(rgN); normalizeWithinSlide(maN, method="l"); rgN <- as.RGData(maN);
rg <- rgN;
append(rg, rgCAvg);
setSlideNames(rg, c("normalized", "raw"));
rm(rgN,rgCAvg);

# Plot M vs A for raw data. From the overlaid grid (according to the
# estimated affine parameters), one can see that the normalization
# method indeed allows some signals to fall out of the positive
# quadrant, that is, the data points outside the grid will become
# negative values (and therefore not visible in the next (A,M)-plot).
ma <- as.MAData(rg);
Alim <- c(0,ceiling(max(ma$A,na.rm=TRUE)));
Mlim <- c(-1,1)*abs(diff(range(Alim)));
plot(ma, slide="raw", xlim=Alim, ylim=Mlim);
drawRGGrid(afit=fit[[1]], drawCurve=TRUE);
plot(ma, slide="normalized", xlim=Alim, ylim=Mlim);
drawRGGrid(maxSignal=Alim[2]);

# Plot M vs M for normalized and raw data. Stratify on observations
# close to A = [0.0,0.5,1.0,...,16.0] to better see the intenstity
# dependency of the log-ratio bias.
ma <- as.MAData(rg);

isCloseTo <- function(x, value=0, dx=0.05) (x > value-dx & x < value+dx);
incl <- rep(NA, nbrOfSpots(rg));
for (aa in seq(0,16, by=0.5))
  incl <- incl | isCloseTo(ma$A[,2], value=aa, dx=0.05);

# Shift the plot to maximize resolution, but keep aspect ratios.
Mlim <- c(-7,4);
plotMvsM(ma, slides=1:2, pch=176, cex=0.5, include=incl, xlim=Mlim);
abline(a=0,b=1, col="lightblue");


}

\seealso{
  \code{\link[aroma.light]{calibrateMultiscan.matrix}}.
  \code{\link{normalizeAffine.RGData}}().
  For more information see \code{\link{RGData}}.
}
\keyword{internal}
\keyword{methods}
