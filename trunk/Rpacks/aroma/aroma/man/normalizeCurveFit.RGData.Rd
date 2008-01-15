%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeCurveFit.RGData}
\alias{normalizeCurveFit.RGData}
\alias{RGData.normalizeCurveFit}
\alias{normalizeCurveFit.RGData}
\alias{normalizeCurveFit,RGData-method}

\alias{normalizeLowess}
\alias{normalizeLoess}
\alias{normalizeSpline}
\alias{normalizeRobustSpline}

\title{Within-slide intensity-dependent normalization in (A,M)}

\usage{\method{normalizeCurveFit}{RGData}(this, slides=NULL, groupBy=NULL, ...)}

\arguments{
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link{LayoutGroups}} specifying the
    groups of spots that are to normalized individually.
    If \code{\link[base]{NULL}}, all spots are normalized together.}
  \item{...}{Other arguments, such as \code{groupBy} and \code{method},
     passed to the curve estimator \code{\link[aroma:estimateMACurve.RGData]{*estimateMACurve}()}.}
}

\value{
  Returns itself invisibly.
}

\description{
  Within-slide intensity-dependent normalization in (A,M). It normalizes pairs of channels by estimating a smooth
  log-ratio intensity-dependent curve.

  Because of they are so well known by their de facto names, the methods
  \code{normalizeLowess(...)} and \code{normalizeLoess(...)} are aliases
  for \code{normalizeCurveFit(..., method="lowess")} and
  \code{normalizeCurveFit(..., method="loess")}, respectively.
}

\section{Negative, non-positive, and saturated values}{
 Non-positive values are set to not-a-number (\code{\link[base:is.finite]{NaN}}).
 Data points that are saturated in one or more channels are not used
 to estimate the normalization function, but they are normalized.
}

\section{Missing values}{
 The estimation of the affine normalization function will be made based
 on only complete non-saturated observations, i.e. observations that
 contains no \code{\link[base]{NA}} values nor saturated values as defined by \code{satSignal}.
}

\section{Weighted normalization}{
 Each data point can be assigned a weight in [0,1] specifying how much
 it should \emph{affect the fitting of the curve-fit normalization function}.
 Note that here a \emph{data point} is here considered to be the pair
 of values in the two channels to be normalized.
 For instance, for two-channel data, a data point is the pair (R,G).

 Regardless of weights, \emph{all} data points are \emph{normalized} based
 on the fitted normalization function.

 Weights can be set using \code{\link[aroma:setProbeWeights.RGData]{*setProbeWeights}()}.
 If weights are specified, they will be used.
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
rg <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
setLayout(rg, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg, c("500V","600V","700V","800V"));
# Reorder the slides in scan order
keepSlides(rg, c("800V","500V","600V","700V"));
# Since scan "500V" seems to be an outlier, exclude it!
removeSlides(rg, slide="500V")

rg0 <- clone(rg);

# Calibrate using x <- (y-a)/b. A better alternative is to projecting the
# fitted line onto (1,1,...,1), but in this example we want to show that
# the scan after translation and rescaling are very similar, but with noise.
fit <- calibrateMultiscan(rg, project=FALSE);

avg <- median
R <- apply(rg$R, MARGIN=1, FUN=avg, na.rm=TRUE)
G <- apply(rg$G, MARGIN=1, FUN=avg, na.rm=TRUE)
R <- R / 8; G <- G/8;
rg <- RGData(R=R,G=G,layout=getLayout(rg))
rm(R,G)

# Since we used project=TRUE all scans are identical now.
keepSlides(rg, slide=1)

Device$set(2, height="108\%");
subplots(9, nrow=3);

# Plot such that there is a right angle between the log(R) and log(G) axes.
Alim <- c(-2,18)
Mlim <- c(-1,1)*abs(diff(Alim))

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main="multiscan calibrated");

# Curve-fit using splines is much faster than using loess and give
# basically the same estimate.
fit <- normalizeCurveFit(rg, method="spline")

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main="calibrated + normalized");
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main="calibrated + normalized");


rg <- clone(rg0)
keepSlides(rg, slide="600V")

ma <- as.MAData(rg);
pmt <- paste("PMT", getSlideNames(ma));
plot(ma, xlim=Alim, ylim=Mlim, main=pmt);

cat("Signals (", pmt, ") before normalization:\n", sep="");
print(summary(rg))
fit <- normalizeCurveFit(rg, method="spline")
cat("Signals (", pmt, ") after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main=paste(pmt, "normalized"));
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main=paste(pmt, "normalized"));


rg <- clone(rg0)

ma <- as.MAData(rg);
plot(NA, xlim=Alim, ylim=Mlim, main="All scans");
for (kk in seq(ma))
  points(ma, slide=kk, col=kk+1, pch=".");

cat("Signals before normalization:\n", sep="");
print(summary(rg))
fit <- normalizeCurveFit(rg, method="spline")
cat("Signals after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)

# To give evidence that all signals now contains the same bias
# and that the relative scale between all scans and all channels
# is one we rescale all signals by the same factor. The only affect
# this will have on the M vs A plot is that the data is shifted
# along the A dimension. Note that there is no rescaling the
# log-ratios!
#scale <- 1/2^5;   # Shift such that A <- A - 5;
#for (ch in c("R", "G"))
#  rg[[ch]] <- scale*rg[[ch]];

ma <- as.MAData(rg);
plot(NA, xlim=Alim, ylim=Mlim, main="All scans normalized");
for (kk in seq(ma))
  points(ma, slide=kk, col=kk+1, pch=".");
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main="All scans normalized");


}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, the light-weight function
  \code{\link[aroma.light]{normalizeCurveFit.matrix}} is used.

  For more information see \code{\link{RGData}}.
}
\keyword{internal}
\keyword{methods}
