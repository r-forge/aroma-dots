fitTpDensity2 <- function(C, M, fpRate=0.01, ..., d=NULL) {
  if (is.null(d)) {
    d <- rep(NA, nrow(M));
  }

  jjs <- whichVector(is.na(d));

  # Done?
  if (length(jjs) == 0)
    return(d);

  # Remove dimnames and translate to make things faster
  dimnames(C) <- NULL;
  dimnames(M) <- NULL;
  C <- t(C);
  M <- t(M);

  for (jj in jjs) {
    if (jj %% 100 == 0)
      print(c(mm, jj, fpRate));
    d[jj] <- .subset2(findTpAtFp(C[,jj], M[,jj], fpRate=fpRate, ..., .checkArgs=FALSE), "tpRate"); 

  } # for (mm in ...)

  # Return updates
  d;
} # fitTpDensity()



findSmoothingForTpAtFp <- function(truth, data, x, fpRate=0.05, minTpRate=0.95, nstepsR=2, accTp=0.001, accR=0.01, ..., verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }
  }

  # Reorder data by position
  o <- order(x);
  x0 <- x[o];

  data0 <- data[o,];
  truth0 <- truth[o,];
  n <- nrow(truth0);

  hTp <- matrix(0, nrow=1, ncol=3);

  h <- c(1,5,10);
  lastR <- Inf;  
  iter <- 1;
  while (diff(range(h)) > accR) {
    verbose && enter(verbose, sprintf("Iteration #%d", iter));
    verbose && printf(verbose, "minTpRate: %.4f\n", minTpRate);
    verbose && printf(verbose, "accTpRate: %.4f\n", accTp);
    verbose && printf(verbose, "Levels: %s\n", paste(h, collapse=","));

    idxs <- getBlockAverageMap(n=n, h=h[2]);
    hApprox2 <- attr(idxs, "hApprox");
    data <- blockAvg(data0, idxs);
    truth <- blockAvg(truth0, idxs);

    fit <- findTpAtFp(truth, data, fpRate=fpRate, ..., verbose=less(verbose), .checkArgs=TRUE);
    tpRate2 <- fit$tpRate;
    verbose && printf(verbose, "tpRate @ %.4f (~%.4f): %.4f\n", h[2], hApprox2, tpRate2);

    if (!h[2] %in% hTp[,1]) {
      t <- c(h[2], hApprox2, tpRate2);
      hTp <- rbind(hTp, t);
      o <- order(hTp[,1]);
      hTp <- hTp[o,];
    }

    if (tpRate2 > minTpRate) {
      # All done?
      if (abs(tpRate2-minTpRate) < accTp)
        break;
      hUse <- h[1];
    } else {
      hUse <- h[3];
    }

    idxs <- getBlockAverageMap(n=n, h=hUse);
    hApprox <- attr(idxs, "hApprox");
    data <- blockAvg(data0, idxs);
    truth <- blockAvg(truth0, idxs);

    fit <- findTpAtFp(truth, data, fpRate=fpRate, ..., verbose=less(verbose), .checkArgs=TRUE);
    tpRate <- fit$tpRate;
    verbose && printf(verbose, "tpRate @ %.4f (~%.4f): %.4f\n", hUse, hApprox, tpRate);

    if (!hUse %in% hTp[,1]) {
      t <- c(hUse, hApprox, tpRate);
      hTp <- rbind(hTp, t);
      o <- order(hTp[,1]);
      hTp <- hTp[o,];
    }

    if (tpRate > minTpRate) {
      # All done?
      if (abs(tpRate-minTpRate) < accTp)
        break;
    }

    if (tpRate2 > minTpRate) {
      h <- c(h[1], h[2]);
      # Approximate estimate; need adjustment 10% of boundary?
      if (tpRate > minTpRate)
        h[1] <- h[1] - 0.05*diff(h);
    } else {
      h <- c(h[2], h[3]);
      # Approximate estimate; need adjustment 10% of boundary?
      if (tpRate < minTpRate)
        h[2] <- h[2] + 0.05*diff(h);
    }
    h <- c(h[1], (h[1]+h[2])/2, h[2]);

    print(hTp);

    iter <- iter + 1;
    verbose && exit(verbose);
  } # while(...)

  hTp;
} # findSmoothingForTpAtFp()



scanTpAtFp <- function(truth, data, x, W=NULL, hs=seq(from=1,to=10,by=0.1), hTp=NULL, ..., shifts=0, adjustForNAs=TRUE, verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }
  }

  verbose && enter(verbose, "Calculating TP rates for different amount of smoothing");

  # Reorder data by position
  o <- order(x);
  x0 <- x[o];

  data0 <- data[o,];
  truth0 <- truth[o,];
  n <- nrow(truth0);

  verbose && printf(verbose, "Data dimensions: %dx%d\n", nrow(data0), ncol(data0));
  if (!identical(dim(data0), dim(truth0)))
    throw("Internal error: ");

  
  nbrOfNAs <- sum(is.na(data0));
  verbose && printf(verbose, "Number of missing values: %d\n", nbrOfNAs);
  
  # Calculate the fraction (0 <= c <= 1) of non-missing values
  fractionOKs <- (1 - nbrOfNAs / length(data0));

  verbose && printf(verbose, "Fraction of OK observations: %.4f%%\n",
                                                     100*fractionOKs);

  if (is.null(hTp)) {
    hTp <- matrix(NA, nrow=length(hs), ncol=5);
    hTp[,1] <- hs;
    colnames(hTp) <- c("h", "hApprox", "xMean", "xMedian", "tpRate");
  } else {
    hs <- hs[!(hs %in% hTp[,1])];
    # Expand hTp for these new ones
    t <- matrix(NA, nrow=length(hs), ncol=5);
    t[,1] <- hs;
    hTp <- rbind(hTp, t);
    rm(t);
    o <- order(hTp[,1]);
    hTp <- hTp[o,];
    rm(o);
  }

  # Skip already existing ones.

  for (kk in seq(along=hs)) {
    h <- hs[kk];

    verbose && enter(verbose, sprintf("Iteration #%d of %d", kk, length(hs)));
    verbose && printf(verbose, "h: %.4f\n", h);
    if (h > 1) {
      data <- truth <- c();
      for (s in shifts) {
        idxs <- getBlockAverageMap(n=n, h=h, s=s);
        if (is.null(W)) {
          data <- rbind(data, blockAvg(data0, idxs));
        } else {
          data <- rbind(data, blockAvg(data0, idxs, 
                              FUN=rowWeightedMeans.matrix, W=W));
        }
        truth <- rbind(truth, blockAvg(truth0, idxs));
      }
#      x <- x / length(shifts);
#      data <- data / length(shifts);
#      truth <- truth / length(shifts);
      if (is.null(W)) {
        x <- blockAvg(x0, idxs);
      } else {
        x <- blockAvg(x0, idxs, FUN=rowWeightedMeans.matrix, W=W);
      }
      hApprox <- attr(idxs, "hApprox");
    } else {
      hApprox <- h;
      x <- x0;
      data <- data0;
      truth <- truth0;
    }
    
    fit <- findTpAtFp(truth, data, ..., verbose=less(verbose), .checkArgs=FALSE);
    tpRate <- fit$tpRate;

    dx <- diff(as.vector(x));

    rr <- whichVector(h == hTp[,1]);

    hh <- c(hApprox, mean(dx, na.rm=TRUE), median(dx, na.rm=TRUE));
    if (adjustForNAs) {
      # Rescale effective amount of smoothing (should become greater...)
      # Rescale distances (should become greater...)
      hh <- hh / fractionOKs;
    }
    hTp[rr,] <- c(h, hh, tpRate);
    verbose && print(verbose, hTp);

    verbose && exit(verbose);
  } # while(...)

  attr(hTp, "fractionOKs") <- fractionOKs;
  attr(hTp, "adjustForNAs") <- adjustForNAs;
  
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  hTp;
} # scanTpAtFp()

############################################################################
# HISTORY:
# 2007-08-20
# o Added file caching to fitRoc2().
# 2007-08-19
# o Renamed argument 'call' to 'toCall' in fitRoc().
# 2007-04-15
# o Added scanTpAtFp().
# 2007-04-14
# o Added interpolation in findTpAtFp().
# o Removed gc() from fitRoc().
# o Added findTpAtFp() to locate TP rate for given FP rate.
# 2007-03-2x
# o Added fitRoc().
# o Created.
############################################################################
