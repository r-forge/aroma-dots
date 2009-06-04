plotDensities <- function(X, xlab=expression(log[2](y[c])), ylab="density", xlim=NULL, ylim=NULL, add=FALSE, ...) {
  if (ncol(X) < 2)
    throw("Argument 'X' must have at least two columns: ", ncol(X));

  nbrOfChannels <- ncol(X);

  dd <- list();
  rangeX <- rangeY <- NA;
  for (ii in seq(length=nbrOfChannels)) {
    x <- X[, ii];
    dd[[ii]] <- density(log(x[is.finite(x) & x > 0], base=2), from=0);
    rangeY <- range(c(rangeY, dd[[ii]]$y, na.rm=TRUE));
    rangeX <- range(c(rangeX, dd[[ii]]$x, na.rm=TRUE));
  }

  if (any(is.na(rangeX)))
    rangeX <- c(0,16);

  if (any(is.na(rangeY)))
    rangeY <- c(0,0.8);

  if (is.null(xlim))
    xlim <- rangeX;

  if (is.null(ylim))
    ylim <- rangeY;

  if (!add)
    plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);

  col <- 1;
  for (ii in seq(length=nbrOfChannels)) {
    line <- dd[[ii]];
    lines(line, col=col, lwd=2, ...);
    col <- col + 1;
  }
} # plotDensities()


