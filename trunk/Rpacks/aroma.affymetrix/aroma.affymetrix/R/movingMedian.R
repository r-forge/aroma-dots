setMethodS3("movingMedian", "matrix", function(Y, x, w=NULL, h=1, na.rm=FALSE, ...) {
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

  # At each position, calculate the weighed average using a 
  # Gaussian kernel.
  for (kk in seq_len(n)) {
    if (kk %% 100 == 0)
      print(kk);

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
}) # movingMedian()



############################################################################
# HISTORY:
# 2007-04-08
# o Added Gaussian smoothing for columns in a matrix.
# 2007-04-02
# o Created.
############################################################################
