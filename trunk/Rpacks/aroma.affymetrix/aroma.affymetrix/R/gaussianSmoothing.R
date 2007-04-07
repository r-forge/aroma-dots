# M_i' = w*M = w*(T-R) = w*T - w*R = T_i' - R'

# Before smoothing, the reference R_i == median(T_i). 
# Keep this property for R' too.

# R' = median(T_i')
# T_i' = M_i' - R'

# => w*T = w*M + w*R = M' + w*R


setMethodS3("gaussianSmoothing", "numeric", function(y, x, w=NULL, sd=1, na.rm=FALSE, ...) {
  # Argument 'y'
  n <- length(y);
  
  # Argument 'x'
  if (length(x) != n) {
    throw("Argument 'y' and 'x' are of different lengths: ", 
                                                 n, " != ", length(x));
  }

  # Argument 'w'
  if (!is.null(w)) {
    if (length(w) != n) {
      throw("Argument 'w' and 'y' are of different lengths: ", 
                                                 length(w), " != ", n);
    }
  }


  # Allocate vector of normalized signals
  yn <- double(n);

  # At each position, calculate the weighed average using a 
  # Gaussian kernel.
  for (kk in 1:n) {
    if (kk %% 100 == 0)
      print(kk);

    # Weights centered around x[kk]
    wKernel <- dnorm(x-x[kk], sd=sd);
    if (!is.null(w))
      wKernel <- w * wKernel;
    wKernel <- wKernel / sum(wKernel);

    # Exclude NAs
    if (na.rm) {
      yn[kk] <- sum(wKernel*y, na.rm=TRUE);
    } else {
      yn[kk] <- sum(wKernel*y);
    }
  }

  yn;
}) # gaussianSmoothing()

############################################################################
# HISTORY:
# 2007-04-02
# o Created.
############################################################################
