setMethodS3("fitProbePositionEffects", "numeric", function(y, seqs, ..., intercept=TRUE, na.rm=TRUE, verbose=FALSE) {
  # Argument 'y':
  y <- Arguments$getDoubles(y, disallow=NULL);

  # Argument 'seqs':
  res <- NULL;
  if (is.character(seqs)) {
    K <- length(seqs);
  } else if (is.raw(seqs)) {
    if (!is.matrix(seqs)) {
      throw("Argument 'seqs' is of type raw, but not a matrix.");
    }
    K <- nrow(seqs);
  } else if (is.list(seqs)) {
    res <- seqs;
    X <- res$X;
    if (!is.matrix(X)) {
      throw("Argument 'seqs' is a list, but does not contain design matrix 'X'.");
    }
    K <- nrow(X);
    rm(X);
    B <- res$B;
    if (!is.matrix(B)) {
      throw("Argument 'seqs' is a list, but does not contain basis matrix 'B'.");
    }
    P <- nrow(B);
    rm(B);

    factors <- res$factors;
    if (!is.vector(factors)) {
      throw("Argument 'seqs' does not contain character vector 'factors'");
    }
  } else {
    throw("Argument 'seqs' is of unknown type: ", class(seqs)[1]);
  }

  if (K != length(y)) {
    throw("Number of probe sequences does not match the number of data points: ", K, " != ", length(y));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && cat(verbose, "Signals:");
  verbose && str(verbose, y);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Excluding missing data points
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (na.rm) {
    # Identify missing sequences
    if (is.character(seqs)) {
      nas <- is.na(seqs);
    } else if (is.raw(seqs)) {
      nas <- (seqs[,1] == as.raw(0));
    } else if (is.list(seqs)) {
      # If passing a design matrix, we assume it does not contain missing
      # values, and if there are, they will add to the estimate of the
      # intercept.
      nas <- FALSE;
    }

    # Identify missing observations
    nas <- nas | is.na(y);

    # Non-missing data points
    nas <- !nas;
    keep <- which(nas);
    rm(nas);

    # Clean out missing data points
    if (length(keep) < K) {
      verbose && enter(verbose, "Exluding missing data points");
      y <- y[keep];
      gc <- gc();

      if (is.character(seqs)) {
        seqs <- seqs[keep];
      } else if (is.raw(seqs)) {
        seqs <- seqs[keep,,drop=FALSE];
      } else if (is.list(seqs)) {
        res$X <- res$X[keep,,drop=FALSE];
      }
      gc <- gc();
    }
    rm(keep);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Building design matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(res)) {
    verbose && enter(verbose, "Building design matrix");

    verbose && printf(verbose, "Sequences (%d):\n", K);
    verbose && str(verbose, seqs);

    res <- getProbePositionEffectDesignMatrix(seqs, ..., 
                            intercept=intercept, verbose=less(verbose, 1));
    verbose && exit(verbose);
  }

  X <- res$X;
  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, X);

  B <- res$B;
  verbose && cat(verbose, "Basis vectors:");
  verbose && str(verbose, B);

  map <- res$map;

  factors <- res$factors;
  verbose && cat(verbose, "Factors:");
  verbose && str(verbose, factors);

  rm(res);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit linear regression model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting model");
  fit <- lm.fit(X, y);
  rm(X, y);
  coefs <- coefficients(fit);
  names(coefs) <- NULL;
  rm(fit);
  gc <- gc();
  verbose && cat(verbose, "Coeffients:");
  verbose && print(verbose, coefs);
  verbose && exit(verbose);

  params <- list();
  if (intercept) {
    params$intercept <- coefs[1];
    coefs <- coefs[-1];
  }

  df <- length(coefs)/length(factors);
  verbose && cat(verbose, "Degrees of freedom: ", df);
  idxs <- 1:df;
  for (kk in seq(along=factors)) {
    key <- names(factors)[kk];
    if (is.null(key)) {
      key <- sprintf("factor%02d", kk);
    }
    params[[key]] <- coefs[idxs];
    coefs <- coefs[-idxs];
  }

  fit <- list(params=params, map=map, B=B);
  class(fit) <- "ProbePositionEffects";

  fit;
}) # fitProbePositionEffects()



setMethodS3("getEffects", "ProbePositionEffects", function(fit, intercept=FALSE, ...) {
  params <- fit$params;
  B <- fit$B;
  map <- fit$map;

  factors <- names(params);
  factors <- setdiff(factors, "intercept");
  F <- length(factors);
  rho <- matrix(0, nrow=nrow(B), ncol=F+1);
  for (kk in 1:F) {
    key <- factors[kk];
    rho[,kk] <- B %*% params[[key]];
  }

  if (intercept) {
    rho <- rho[,1:F];
    colnames(rho) <- factors;
  } else {
    rho <- rho - rowSums(rho[,1:F])/(F+1);
    colnames(rho) <- names(map)[-1];
  }

  rho;
}) # getEffects()



setMethodS3("predict", "ProbePositionEffects", function(object, seqs, ..., verbose=FALSE) {
  # To please R CMD check
  fit <- object;


  # Argument 'seqs':
  if (is.character(seqs)) {
    K <- length(seqs);
    P <- nchar(seqs);
    P <- unique(P);
    if (length(P) != 1) {
      throw("Argument 'seqs' contains sequences of different lengths: ",
                            paste(head(sort(P)), collapse=", "));
    }
  } else if (is.raw(seqs)) {
    if (!is.matrix(seqs)) {
      throw("Argument 'seqs' is of type raw, but not a matrix.");
    }
    K <- nrow(seqs);
  } else {
    throw("Argument 'seqs' is of unknown type: ", class(seqs)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Predicting probe-affinities from probe-position parameters");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coercing sequences into a raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(seqs)) {
    verbose && enter(verbose, "Coercing sequences into a raw sequence matrix");
    K <- length(seqs);
    P <- nchar(seqs[1]);
    verbose && cat(verbose, "Number of sequences: ", K);
    verbose && cat(verbose, "object.size(seqs):");
    verbose && print(verbose, object.size(seqs));
    seqs <- paste(seqs, collapse="");
    seqs <- charToRaw(seqs);
    throw("Not implemented yet, i.e. to be fixed.");
    seqs <- matrix(seqs, nrow=K, ncol=P, byrow=TRUE);
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && cat(verbose, "object.size(seqs):");
    verbose && print(verbose, object.size(seqs));
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get probe-position effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rho <- getEffects(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate predicted values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  K <- nrow(seqs);
  P <- ncol(seqs);

  # Allocate probe-affinity vector
  phi <- double(K);

  map <- attr(seqs, "map");
  values <- map[-1];
  factors <- names(values);
#  values <- map[c("A", "C", "G", "T")];

  # For each position
  for (pp in seq(length=P)) {
    verbose && enter(verbose, sprintf("Probe position #%d of %d", pp, P));

    # Get the nucleotides at this position for all sequences
    seqsPP <- seqs[,pp];

    for (bb in 1:ncol(rho)) {
#      verbose && enter(verbose, sprintf("Factor #%d ('%s') of %d", bb, factors[bb], ncol(rho)));

      # Identify sequences with nucleotide 'bb' at position 'pp'.
      idxs <- which(seqsPP == values[bb]);

      # Add the nucleotide effect rho(pp,bb) to the probe-affinity
      phi[idxs] <- phi[idxs] + rho[pp,bb];

#      verbose && exit(verbose);
    } # for (bb ...)

    rm(seqsPP);
    verbose && exit(verbose);
  } # for (pp ...)

  verbose && exit(verbose);

  phi;
}) # predict()



setMethodS3("plot", "ProbePositionEffects", function(x, type="b", col=NULL, pch=NULL, lwd=2, xlab="Position", ylab="Effect (on log2 scale)", ...) {
  # To please R CMD check
  fit <- x;

  rho <- getEffects(fit);
  if (is.null(col)) {
    col <- seq(length=ncol(rho));
  }
  if (is.null(pch)) {
    pch <- colnames(rho);
    if (any(nchar(pch) > 1)) {
      pch <- seq(along=pch);
    }
  }
  matplot(rho, type=type, 
          col=col, lwd=lwd, pch=pch, 
          xlab=xlab, ylab=ylab, ...);
  abline(h=0, lty=3);
  
#  legend("topleft", colnames(rho), pch=pch, col=col);
}) # plot();




setMethodS3("text", "ProbePositionEffects", function(x, labels=NULL, col=NULL, ...) {
  # To please R CMD check
  fit <- x;

  rho <- getEffects(fit);

  if (is.null(labels)) {
    labels <- colnames(rho);
  }

  if (is.null(col)) {
    col <- seq(length=ncol(rho));
  }

  xx <- seq(length=nrow(rho));
  for (cc in seq(length=ncol(rho))) {
    yy <- rho[,cc];
    text(xx,yy, labels=labels[cc], col=col[cc], ...);
  }
}) # plot();






############################################################################
# HISTORY:
# 2008-07-11
# o Renamed getProbePositionEffects() to getEffects().
# o Updated to work with the new AromaCellSequenceFile encodings.
# 2008-07-07
# o Added predict().
# o Added getProbePositionEffects().
# o Added na.rm=TRUE to fitProbePositionEffects().
# 2008-07-05
# o Created.
############################################################################
