%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeFragmentLength.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{normalizeFragmentLength}
\alias{normalizeFragmentLength.default}
\alias{normalizeFragmentLength}

 \title{Normalizes signals for PCR fragment-length effects}

 \description{
  Normalizes signals for PCR fragment-length effects. Some or all signals are used to estimated the
  normalization function.  All signals are normalized.
 }

 \usage{\method{normalizeFragmentLength}{default}(y, fragmentLengths, targetFcns=NULL, subsetToFit=NULL, .isLogged=TRUE, ..., .returnFit=FALSE)}

 \arguments{
   \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length K of signals to be normalized
     across E enzymes.}
   \item{fragmentLengths}{An \code{\link[base]{integer}} KxE \code{\link[base]{matrix}} of fragment lengths.}
   \item{targetFcns}{A \code{\link[base]{list}} of E \code{\link[base]{function}}s - one per enzyme.}
   \item{subsetToFit}{The subset of data points used to fit the
      normalization function.
      If \code{\link[base]{NULL}}, all data points are considered.}
   \item{.isLogged}{A \code{\link[base]{logical}}.}
   \item{...}{Additional arguments passed to \code{\link[stats]{lowess}}.}
   \item{.returnFit}{A \code{\link[base]{logical}}.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of the normalized signals.
 }

 \section{Multi-enzyme normalization}{
  It is assumed that the fragment-length effects from multiple enzymes
  added (with equal weights) on the intensity scale.
  The fragment-length effects are fitted for each enzyme separately based
  on units that are exclusively for that enzyme.
  \emph{If there are no or very such units for an enzyme, the assumptions
  of the model are not met and the fit will fail with an error.}
  Then, from the above single-enzyme fits the average effect across
  enzymes is the calculated for each unit that is on multiple enzymes.
 }

 \examples{
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 1: Single-enzyme fragment-length normalization of 6 arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number samples
I <- 9;

# Number of loci
J <- 1000;

# Fragment lengths
fl <- seq(from=100, to=1000, length.out=J);

# Simulate data
y <- matrix(0, nrow=J, ncol=I);
maxY <- 12;
for (kk in 1:I) {
  k <- runif(n=1, min=3, max=5);
  mu <- function(fl) {
    maxY - fl^{1/k};
  }
  eps <- rnorm(J, mean=0, sd=1);
  y[,kk] <- mu(fl) + eps;
}

# Normalize data (to a zero baseline)
yN <- apply(y, MARGIN=2, FUN=function(y) {
  normalizeFragmentLength(y, fragmentLengths=fl);
})

# Plot raw data
layout(matrix(1:9, ncol=3));
xlim <- c(0,max(fl));
ylim <- c(0,max(y));
xlab <- "Fragment length";
ylab <- expression(log2(theta));
for (kk in 1:I) {
  plot(fl, y[,kk], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
  lines(lowess(fl, y[,kk]), col="red", lwd=2);
}

# Plot normalized data
layout(matrix(1:9, ncol=3));
ylim <- c(-1,1)*max(y)/2;
for (kk in 1:I) {
  plot(fl, yN[,kk], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
  lines(lowess(fl, yN[,kk]), col="blue", lwd=2);
}


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 2: Two-enzyme fragment-length normalization of 6 arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number samples
I <- 3;

# Number of loci
J <- 1000;

# Fragment lengths (two enzymes)
fl <- matrix(0, nrow=J, ncol=2);
fl[,1] <- seq(from=100, to=1000, length.out=J);
fl[,2] <- seq(from=1000, to=100, length.out=J);

# Let only 1/3 of the units be on both enzymes
fl[seq(from=1, to=J, by=3),1] <- NA;
fl[seq(from=2, to=J, by=3),2] <- NA;

# Simulate data
z <- array(0, dim=c(J,2,I));
maxY <- 12;
for (ii in 1:I) {
  # Common effect for both enzymes
  mu <- function(fl) {
    k <- runif(n=1, min=3, max=5);
    maxY - fl^{1/k};
  }
  for (ee in 1:2) {
    eps <- rnorm(J, mean=0, sd=1/sqrt(2));
    z[,ee,ii] <- mu(fl[,ee]) + eps;
  }
}
# The observed data is a mix of two enzymes
y <- matrix(NA, nrow=J, ncol=I);
hasFl <- is.finite(fl);
# Single-enzyme units
ok <- (hasFl[,1] & !hasFl[,2]);
y[ok,] <- z[ok,1,];
ok <- (!hasFl[,1] & hasFl[,2]);
y[ok,] <- z[ok,2,];
# Both-enzyme units (sum on intensity scale)
ok <- (hasFl[,1] & hasFl[,2]);
y[ok,] <- log2((2^z[ok,1,]+2^z[ok,2,])/2);

# Normalize data (to a zero baseline)
yN <- matrix(NA, nrow=J, ncol=I);
fits <- vector("list", I);
for (ii in 1:I) {
  yNi <- normalizeFragmentLength(y[,ii], fragmentLengths=fl, .returnFit=TRUE);
  fits[[ii]] <- attr(yNi, "modelFit");
  yN[,ii] <- yNi;
}

# Plot raw data
layout(matrix(1:(2*I), ncol=3, byrow=TRUE));
xlim <- c(0,max(fl,na.rm=TRUE));
ylim <- c(0,max(y,na.rm=TRUE));
xlab <- "Fragment length";
ylab <- expression(log2(theta));
for (ii in 1:I) {
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main="raw");

  # Single-enzyme units
  for (ee in 1:2) {
    # The raw data
    points(fl[,ee], z[,ee,ii], col=ee+1);
  }

  # Both-enzyme units (use fragment-length for enzyme #1)
  ok <- is.finite(fl[,1]) & is.finite(fl[,2]);
  points(fl[ok,1], y[ok,ii], col=3+1);

  for (ee in 1:2) {
    # The true effects
    ok <- is.finite(fl[,ee]);
    lines(lowess(fl[ok,ee], z[ok,ee,ii]), col="black", lwd=3);

    # The estimated effects
    fit <- fits[[ii]][[ee]]$fit;
    lines(fit, col="orange", lwd=2);
  }
}

# Plot normalized data
ylim <- c(-1,1)*max(y,na.rm=TRUE)/2;
ok <- !(hasFl[,1] & hasFl[,2]);
for (ii in 1:I) {
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main="normalized");
  # Single-enzyme units
  for (ee in 1:2) {
    # The normalized data
    ok <- is.finite(fl[,ee]) & !is.finite(fl[,-ee]);
    points(fl[ok,ee], yN[ok,ii], col=ee+1);
  }
  # Both-enzyme units (use fragment-length for enzyme #1)
  ok <- is.finite(fl[,1]) & is.finite(fl[,2]);
  points(fl[ok,1], yN[ok,ii], col=3+1);
}


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
