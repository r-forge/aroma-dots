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

