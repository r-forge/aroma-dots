# rowMedians.matrix() for now, because it support weights(!)
# (although 'matrixStats' exists).
library(R.native);

getBlockAverageMap <- function(n, h=1, s=0, ...) {
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(1,1000));


  # Is h an integer?
  if (h == as.integer(h)) {
    idxs <- matrix(seq_len(n), nrow=h);
    if (n %% h != 0)
      idxs <- idxs[,-ncol(idxs)];
  } else {
    K <- ceiling(n/h);
    idxs <- matrix(TRUE, nrow=ceiling(h), ncol=K);
    steps <- (h %% 1) * ceiling(K);
    incl <- seq(from=1, to=K, length=steps);
    incl <- round(incl);
    idxs[ceiling(h), -incl] <- FALSE;
    # Shift
    if (s > 0) {
      lastRow <- idxs[ceiling(h),];
      tail <- seq(from=length(lastRow)-s+1, to=length(lastRow));
      lastRow <- c(lastRow[tail], lastRow[-tail]);
      idxs[ceiling(h),] <- lastRow;
    }
    idxs[idxs] <- seq_len(n);
    idxs[idxs == 0] <- NA;
  }

  # Skip last column in case looping
  if (n %% h != 0)
    idxs <- idxs[,-ncol(idxs)];

  # The effective 'h'
  hApprox <- sum(!is.na(idxs))/ncol(idxs);
  attr(idxs, "hApprox") <- hApprox;

  idxs;
} # getBlockAverageMap()


blockAvg <- function(X, idxs, FUN=rowMedians.matrix, W=NULL, ...) {
  na.rm <- (any(is.na(X)) || any(is.na(idxs)));
  dimnames <- dimnames(X);
  dimnames(X) <- NULL;
  X <- t(X);
  if (!is.null(W))
    W <- t(W);

  Z <- apply(idxs, MARGIN=2, FUN=function(jj) {
    jj <- jj[is.finite(jj)];
    Zjj <- X[,jj,drop=FALSE];
    if (!is.null(W)) {
      Wjj <- W[,jj,drop=FALSE];
      Zjj <- FUN(Zjj, W=Wjj, ..., na.rm=na.rm);
    } else {
      Zjj <- FUN(Zjj, W=W, ..., na.rm=na.rm);
    }
    Zjj;        
  });
  Z <- t(Z);
  colnames(Z) <- dimnames[[2]];
  Z;
} # blockAvg()


##############################################################################
# HISTORY:
# 2008-01-11
# o Extracted from CRMA-Fig8,res,filtered.R.
# 2007-11-16
# o Added argument 'W' to blockAvg().
############################################################################## 
