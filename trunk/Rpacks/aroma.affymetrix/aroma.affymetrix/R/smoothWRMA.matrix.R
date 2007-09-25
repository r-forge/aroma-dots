setMethodS3("smoothWRMA", "matrix", function(Y, x, kernel=gaussKernel, sd=100e3, na.rm=TRUE, ..., progress=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'Y'
  K <- nrow(Y);  # Number of positions
  I <- ncol(Y);  # Number of samples
  
  # Argument 'x'
  if (length(x) != K) {
    throw("Argument 'x' has different number of values that rows in 'Y': ", 
                                                     length(x), " != ", K);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit everything on the log2 scale
  Y <- log2(Y);

  # Identify missing values?
  if (na.rm) {
    nas <- is.na(Y);
    dim(nas) <- c(K,I);
  }

  # Allocate vector of smoothed signals
  theta <- matrix(NA, nrow=K, ncol=I);
  phi <- rep(NA, K);

  # At each position, calculate the weighed average using a 
  # Gaussian kernel.
  for (kk in seq(length=K)) {
    if (progress && kk %% 100 == 0)
      print(kk);

    # Weights centered around x[kk]
    w <- kernel(x, mean=x[kk], sd=sd);
    # Weight matrix
    w <- matrix(w, nrow=K, ncol=I);
   
    # Give missing values zero weight?
    if (na.rm)
      w[nas] <- 0;

    wR <- rowSums(w);
    keep <- which(wR > 0);
    if (length(keep) > 0) {
      w <- w[keep,,drop=FALSE];
      y <- Y[keep,,drop=FALSE];
      verbose && print(verbose, list(y=y, w=w));

      # Fit weighted RMA (don't have to rescale weights)
      # Call fitWRMA.matrix() explicitly to avoid dispatching
      fit <- fitWRMA.matrix(y=y, w=w, .log2=FALSE);
      verbose && print(verbose, fit);
      theta[kk,] <- .subset2(fit, "theta");
      # Average phi?
      phi[kk] <- .subset2(fit, "avgPhi");
    } else {
      # If no data, keep theta:s and phi:s as missing values.
    }
  }

  # Return everything on the intensity scale
  theta <- 2^theta;
  phi <- 2^phi;

  list(theta=theta, phi=phi);
}) # smoothWRMA()


############################################################################
# HISTORY:
# 2007-09-24
# o Now smoothWRMA() returns (theta, phi) on the intensity scale.
# 2007-09-18
# o Created from gaussianSmoothing.R.
############################################################################
