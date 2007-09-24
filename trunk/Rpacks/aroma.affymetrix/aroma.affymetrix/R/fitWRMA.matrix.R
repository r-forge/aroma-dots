setMethodS3("fitWRMA", "matrix", function(y, w, psiCode=0, psiK=1.345) {
  I <- ncol(y);
  K <- nrow(y);

  fit <- .Call("R_wrlm_rma_default_model", y, psiCode, psiK, w, PACKAGE="affyPLM");

  est <- fit$Estimates;
  se <- fit$StdErrors;

  # Chip effects
  beta <- est[1:I];

  # Probe affinities
  alpha <- est[(I+1):length(est)];
  alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)]);

  # Estimates on the intensity scale
  theta <- 2^beta;
  phi <- 2^alpha;

  # Calcuate log-ratios
  thetaR <- median(theta, na.rm=TRUE);
  M <- log2(theta/thetaR);

  list(theta=theta, phi=phi, M=M);
}, protected=TRUE) # fitWRMA()


############################################################################
# HISTORY:
# 2007-09-18
# o Created.
############################################################################
