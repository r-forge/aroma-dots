library(ROC);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fitRoc <- function(truth, data, ncuts=2000, cutpts=NULL, call=1, ..., .checkArgs=TRUE) {
  # Call females?
  if (call == 1) {
    # Default rule: function(x, thresh) ifelse(x > thresh, 1, 0);
    # That is, call female in case *above* threshold
  } else {
    # Complement rule: function(x, thresh) ifelse(x < thresh, 1, 0);
    # That is, call male in case *below* threshold
    truth <- 1-truth;
    data <- -data;
  }

  if (.checkArgs) {
    truth <- as.vector(truth);
    data <- as.vector(data);
    if (length(truth) != length(data)) {
      throw("Argument 'truth' and 'data' are of different lengths: ");
    }
  
    ok <- (is.finite(truth) & is.finite(data));
    truth <- truth[ok];
    data <- data[ok];
    rm(ok);
  }

  if (is.null(cutpts)) {
    probs <- seq(from=0, to=1, length=ncuts);
    cutpts <- quantile(data, probs=probs, na.rm=TRUE);
  }

  fit <- rocdemo.sca(truth=truth, data=data, cutpts=cutpts, ...);

  fit@markerLabel <- xlab;
  fit@caseLabel <- ylab;
  if (call == 0) {
    fit@cuts <- -fit@cuts;
  }

  fit;
} # fitRoc()



findTpAtFp <- function(truth, data, fpRate=0.05, nsteps=50, acc=1e-5, ..., verbose=FALSE, .checkArgs=TRUE) {
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


  # Vectorize
  data <- as.vector(data);
  truth <- as.vector(truth);

  ok <- (is.finite(truth) & is.finite(data));
  n <- length(ok);
  nok <- sum(ok);
  if (nok < n) {
    data <- data[ok];
    truth <- truth[ok];
  }
  bad <- which(!ok);

  res <- list(
    tpRate=NA,
    fpRate=NA,
    threshold=NA,
    acc=acc,
    tpRateRange=NA,
    fpRateRange=NA,
    thresholdRange=NA,
    excludedNAs=bad
  );

  if (nok == 0) {
    return(res);
  }

  verbose && enter(verbose, "Finding TP rate at FP rate...");  
  r <- range(data, na.rm=TRUE);

  lastR <- c(-Inf,+Inf);
  while (all(!is.na(r)) && (diff(r) > 2*acc) && !identical(r, lastR)) {  
    verbose && printf(verbose, "FP range: [%.2g,%.2g]\n", r[1], r[2]);  
    cutpts <- seq(from=r[1], to=r[2], length=nsteps);
    fit <- fitRoc(truth=truth, data=data, cutpts=cutpts, ..., .checkArgs=FALSE);
    betas <- 1-fit@spec;
    isLess <- which(fpRate <= betas);
    from <- isLess[length(isLess)];
    isMore <- which(betas <= fpRate);
    to <- isMore[1];
    lastR <- r;
    r <- cutpts[c(from,to)];
    r <- c(r[1], r[2]);  # In case length(r) == 1
    verbose && printf(verbose, "Next FP range: [%.2g,%.2g]\n", r[1], r[2]);  
  }
  verbose && exit(verbose);

  fpRateRange <- range(1-fit@spec, na.rm=TRUE);
  tpRateRange <- range(fit@sens, na.rm=TRUE);
  res$tpRateRange <- tpRateRange;
  res$fpRateRange <- fpRateRange;
  res$thresholdRange <- r;

  if (any(is.na(r))) {
    # Case a) TP rate is zero
    res$tpRate <- 0;
    res$fpRate <- fpRate;
    res$threshold <- mean(r, na.rm=TRUE);
  } else if (fpRateRange[1] == fpRateRange[2]) {
    # Case b) Exact (FP,TP) rate
    res$rho <- 0;
    res$fpRate <- fpRateRange[1];
    res$tpRate <- fit@sens[1];
    res$threshold <- r[1];
  } else {
    # Case c) Interpolate (FP,TP) rate
    rho <- (fpRate-fpRateRange[1])/(fpRateRange[2]-fpRateRange[1]);
    res$rho <- rho;
    w2 <- rho;
    w1 <- 1-w2;
    res$fpRate <- w1*fpRateRange[1] + w2*fpRateRange[2]; 
    res$tpRate <- w1*tpRateRange[1] + w2*tpRateRange[2];
    res$threshold <- r[1] + rho*(r[2]-r[1]);
  }


  res;
} # findTpAtFp()


fitTpDensity <- function(C, M, fpRate=0.01, ..., d=NULL) {
  if (is.null(d)) {
    d <- rep(NA, nrow(M));
  }

  jjs <- which(is.na(d));

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

############################################################################
# HISTORY:
# 2007-04-14
# o Added interpolation in findTpAtFp().
# o Removed gc() from fitRoc().
# o Added findTpAtFp() to locate TP rate for given FP rate.
# 2007-03-2x
# o Added fitRoc().
# o Created.
############################################################################
