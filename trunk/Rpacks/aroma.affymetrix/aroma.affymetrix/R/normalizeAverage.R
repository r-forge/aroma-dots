setMethodS3("normalizeAverage", "list", function(x, baseline=1, avg=median, targetAvg=2200, ...) {
  # Estimate the scale for each channel
  scale <- lapply(x, FUN=avg, ...);
  scale <- unlist(scale, use.names=FALSE);
  scale1 <- scale[baseline];

  # Standardize so that the 'baseline' channel has scale one.
  scale <- scale / scale1;

  # Rescale to target averages?
  if (!is.null(targetAvg)) {
    rho <- (scale1 / targetAvg);
    scale <- rho * scale;
  }

  # Rescale so that all channels have the same scale
  for (cc in 1:length(x)) {
    x[[cc]] <- x[[cc]] / scale[[cc]];
  }

  x;
}, private=TRUE)


setMethodS3("normalizeAverage", "matrix", function(x, baseline=1, avg=median, targetAvg=2200, ...) {
  # Estimate the scale for each channel
  scale <- apply(x, MARGIN=2, FUN=avg, ...);
  scale1 <- scale[baseline];

  # Standardize so that the 'baseline' channel has scale one.
  scale <- scale / scale1;

  # Rescale to target averages?
  if (!is.null(targetAvg)) {
    rho <- (scale1 / targetAvg);
    scale <- rho * scale;
  }

  # Rescale so that all channels have the same scale
  for (cc in 1:ncol(x)) {
    x[,cc] <- x[,cc] / scale[cc];
  }

  x;
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-05-08
# o Created.
############################################################################
