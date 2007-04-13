# M_i' = w*M = w*(T-R) = w*T - w*R = T_i' - R'

# Before smoothing, the reference R_i == median(T_i). 
# Keep this property for R' too.

# R' = median(T_i')
# T_i' = M_i' - R'

# => w*T = w*M + w*R = M' + w*R


setMethodS3("gaussianSmoothing", "matrix", function(Y, x, w=NULL, sd=1, na.rm=FALSE, ...) {
  # Argument 'Y'
  n <- nrow(Y);
  k <- ncol(Y);
  
  # Argument 'x'
  if (length(x) != n) {
    throw("Argument 'x' has different number of values that rows in 'Y': ", 
                                                     length(x), " != ", n);
  }

  # Argument 'w'
  if (!is.null(w)) {
    if (length(w) != n) {
      throw("Argument 'w' has different number of values that rows in 'Y': ", 
                                                       length(w), " != ", n);
    }
  }


  # Allocate vector of smoothed signals
  Ys <- matrix(0, nrow=n, ncol=k);

#  wKernelMax <- dnorm(x, sd=sd);

  # At each position, calculate the weighed average using a 
  # Gaussian kernel.
  for (kk in seq_len(n)) {
    if (kk %% 100 == 0)
      print(kk);

    # Weights centered around x[kk]
    wKernel <- dnorm(x, mean=x[kk], sd=sd);
    wKernel <- wKernel / sum(wKernel);

    # Exclude NAs
    if (na.rm) {
      Ys[kk,] <- colSums(wKernel*Y, na.rm=TRUE);
    } else {
      Ys[kk,] <- colSums(wKernel*Y);
    }
  }

  if (!is.null(w)) {
    w <- w / sum(w, na.rm=TRUE);
    Ys <- w*Ys;
  }

  Ys;
}) # gaussianSmoothing()




setMethodS3("gaussianSmoothing", "numeric", function(y, ...) {
  as.vector(gaussianSmoothing(as.matrix(y), ...));
}) # gaussianSmoothing()


############################################################################
# HISTORY:
# 2007-04-08
# o Added Gaussian smoothing for columns in a matrix.
# 2007-04-02
# o Created.
############################################################################
