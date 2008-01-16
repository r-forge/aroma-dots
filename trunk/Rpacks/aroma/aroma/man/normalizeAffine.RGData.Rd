%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.NORM.R
% 
% on Tue Jan 15 18:36:04 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeAffine.RGData}
\alias{normalizeAffine.RGData}
\alias{RGData.normalizeAffine}
\alias{normalizeAffine.RGData}
\alias{normalizeAffine,RGData-method}


 \encoding{latin1}

 \title{Weighted affine normalization between channels and arrays}

 \description{
   Weighted affine normalization between channels and arrays.

   This method will both remove curvature in the M vs A plots that are
   due to an affine transformation of the data. In other words, if there
   are (small or large) biases in the different (red or green) channels,
   biases that can be equal too, you will get curvature in the M vs A plots
   and this type of curvature will be removed by this normalization method.

   Moreover, if you normalize all slides at once (recommended), this
   method will also bring the signals on the same scale such that the
   log-ratios for different slides are comparable. Thus, do not normalize
   the scale of the log-ratios between slides afterward.

   It is recommended to normalize as many slides as possible in one run.
   The result is that if creating log-ratios between any channels and any
   slides, they will contain as little curvature as possible.

   Furthermore, since the relative scale between any two channels on any
   two slides will be one if one normalizes all slides (and channels) at
   once it is possible to add or multiply with the \emph{same} constant
   to all channels/arrays without introducing curvature. Thus, it is
   easy to rescale the data afterwards as demonstrated in the example.
 }

 \usage{\method{normalizeAffine}{RGData}(this, slides=NULL, groupBy=NULL, ...)}

 \arguments{
  \item{slides}{\code{\link[base]{vector}} of slides to be normalized \emph{at once}.
   If \code{\link[base]{NULL}}, all slides are used.
   To normalize several slides but slide by slide, call this method
   once for each slide instead.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link[aroma]{LayoutGroups}}
   specifying the groups of spots that are to normalized individually.
   If \code{\link[base]{NULL}}, all spots are considered.}
  \item{...}{Additional arguments accepted by
   \code{\link[aroma.light]{normalizeAffine.matrix}}.}
 }

 \value{
   Returns a \code{\link[base]{list}} containing estimated biases, slopes, and distances
   to the optimal line for each slide.
 }

 \section{Negative, non-positive, and saturated values}{
  Affine normalization applies equally well to negative values. Thus,
  contrary to normalization methods applied to log-ratios, such as curve-fit
  normalization methods, affine normalization, will not set these to \code{\link[base]{NA}}.

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
  it should \emph{affect the fitting of the affine normalization function}.
  Note that here a \emph{data point} is here considered to be the \code{\link[base]{vector}}
  of values for all channels and all arrays included in the normalization.
  For instance, for two-channel data with three arrays, each data point is
  a vector (R1,G1,R2,G2,R3,G3) of length 6 where R and G represent the
  red and the green channels.

  Regardless of weights, \emph{all} data points are
  \emph{normalized} based on the fitted normalization function.

  Data-point weights are obtained from probe weights, if given.
  Weights can be set using \code{\link[aroma:setProbeWeights.RGData]{*setProbeWeights}()}.
  If weights are specified, they will be used.
 }

 \section{Within-slide normalization}{
   This normalization method normalized multiple channels/arrays at once.
   To normalize array by array, like curve-fit normalization, use a loop,
   e.g. \code{for (kk in 1:nbrOfSlides(rg)) normalizeAffine(rg, slide=kk)}.
 }

 \details{
  A line is fitted robustly throught the \eqn{(y_R,y_G)} observations
  using an iterated re-weighted principal component analysis (IWPCA),
  which minimized the residuals that are orthogonal to the fitted line.
  Each observation is down-weighted by the inverse of the absolute
  residuals, i.e. in \eqn{L_1}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   [1] Henrik Bengtsson and Ola Hössjer, \emph{Methodological Study of Affine Transformations of Gene Expression Data}, Methodological study of affine transformations of gene expression data with proposed robust non-parametric multi-dimensional normalization method, BMC Bioinformatics, 2006, 7:100.
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
# Reorder the slides in scan order
keepSlides(rg0, c("800V","500V","600V","700V"));
# Since scan "500V" seems to be an outlier, exclude it!
removeSlides(rg0, slide="500V")

rg <- clone(rg0);

# Multiscan calibration.
fit <- calibrateMultiscan(rg);


Device$set(3, height="108\%");
subplots(9, nrow=3);

# Plot such that there is a right angle between the log(R) and log(G) axes.
Alim <- c(-2,18)
Mlim <- c(-1,1)*abs(diff(Alim))

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main="multiscan calibrated");

fit <- normalizeAffine(rg)
drawRGGrid(maxSignal=Alim[2], afit=fit[[1]])

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
fit <- normalizeAffine(rg)
cat("Signals (", pmt, ") after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)
drawRGGrid(maxSignal=Alim[2], afit=fit[[1]])

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
fit <- normalizeAffine(rg)
cat("Signals after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)
drawRGGrid(maxSignal=Alim[2], afit=fit[[1]])

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

 \seealso{
   Internally, the light-weight function
   \code{\link[aroma.light]{normalizeAffine.matrix}} is used.

   For more information see \code{\link{RGData}}.
 }
\keyword{internal}
\keyword{methods}
