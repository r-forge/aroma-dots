###########################################################################/**
# @RdocClass MatrixRocData
#
# @title "The MatrixRocData class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{data}{A @numerical JxK @matrix.}
#   \item{truth}{A @list of K @functions, a @numerical JxK @matrix, or 
#     a @vector of either length J or length K.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("MatrixRocData", function(data=NULL, truth=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'truth':
  if (inherits(truth, "MatrixRocData")) {
    object <- truth;
    truth <- object$.truth;
    data <- object$.data;
  }

  if (!is.null(data)) {
    # Argument 'data':
    if (!is.numeric(data)) {
      throw("Argument 'data' is not numeric: ", mode(data));
    }
    if (!is.matrix(data)) {
      data <- as.matrix(data);
    }
    dim <- dim(data);
    nbrOfEntries <- dim[1];
    nbrOfChannels <- dim[2];

    if (nbrOfChannels == 0) {
      throw("Argument 'data' contains no channels.");
    }

    if (nbrOfEntries < 2) {
      throw("Argument 'data' contains less than two entries: ", nbrOfEntries);
    }


    # Argument 'truth':
    if (is.function(truth)) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # The truth is specified as a global function
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    } else if (is.list(truth)) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # The truth is specified as one function per channel
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (length(truth) != nbrOfChannels) {
        throw("The number of functions in argument 'truth' (list) does not match the number of channels: ", length(truth), " != ", nbrOfChannels);
      }
      for (kk in seq(along=truth)) {
        if (!is.function(truth[[kk]])) {
          throw("Element #", kk, " of argument 'truth' (list) is not a function: ", mode(truth[[kk]]));
        }
      }
    } else if (is.numeric(truth)) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # The truth is specified explicitly and mapped to data points
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # A K vector specified?
      dimT <- dim(truth);
      if (is.null(dimT)) {
        # Infer whether the truth is by entry or by channel
        if (length(truth) == nbrOfEntries) {
          dimT <- c(nbrOfEntries, 1);
        } else if (length(truth) == nbrOfChannels) {
          dimT <- c(1, nbrOfChannels);
        } else {
          throw("The length of argument 'truth' (vector) matches neither the number of entries (rows) nor the number of channels (columns) of 'data': ", length(truth), " not in (", nbrOfEntries, ", ", nbrOfChannels, ")");
        }
        dim(truth) <- dimT;
      }
  
      truth <- as.integer(truth);
      dim(truth) <- dimT;
  
  
      # Expand by entries or channels?
      if (any(dimT != dim)) {
        # By entries or by channels?
        byrow <- (dimT[2] == dim[2]);
  
        # Assert compatibility of dimensions
        cc <- ifelse(byrow, 1, 2);
        if (dimT[cc] != dim[cc]) {
         throw("The dimension of argument 'truth' is not compatible with that of argument 'data': ", paste(dimT, collapse="x"), " != ", paste(dim, collapse="x"));
        }
  
        # Expand
        truth <- matrix(truth, nrow=dim[1], ncol=dim[2], byrow=byrow);
      }
    } else {
      throw("Argument 'truth' must be either a list or a numeric: ", 
                                                         mode(truth));
    }
  }

  extend(RocData(truth=truth, data=data, ...), "MatrixRocData");
}) # MatrixRocData()


setMethodS3("getNEUTRAL", "MatrixRocData", function(this, ...) {
  as.integer(0);
})

setMethodS3("getLOSS", "MatrixRocData", function(this, ...) {
  as.integer(-1);
})

setMethodS3("getGAIN", "MatrixRocData", function(this, ...) {
  as.integer(1);
})

setMethodS3("getUNKNOWN", "MatrixRocData", function(this, ...) {
  as.integer(NA);
})


setMethodS3("getData", "MatrixRocData", function(this, ...) {
  Y <- this$.data;
  Y;
})

setMethodS3("setData", "MatrixRocData", function(this, Y, ...) {
  oldY <- getData(this);
  # Argument 'Y':
  if (any(dim(Y) != dim(oldY))) {
    throw("Cannot set new data. Incompatible dimensions: ", paste(dim(Y), collapse="x"), " != ", paste(dim(oldY), collapse="x"));
  }

  this$.data <- Y;

  invisible(this);
})


setMethodS3("getTruth", "MatrixRocData", function(this, byNames=TRUE, ...) {
  T <- this$.truth;

  # Apply functions?
  if (is.function(T) || is.list(T)) {
    x <- getPositions(this);
    Y <- getData(this);

    if (is.list(T)) {
      fcns <- T;
    } else {
      fcns <- rep(list(T), ncol(Y));
    }

#    channels <- colnames(Y);
#    if (is.null(channels))
    channels <- seq(length=ncol(Y));

    # Allocate results
    T <- matrix(0L, nrow=nrow(Y), ncol=ncol(Y));
    for (kk in seq(length=ncol(Y))) {
      fcn <- fcns[[kk]];
      t <- fcn(x, channel=channels[kk], y=Y[,kk], ...);
      if (is.character(t)) {
        t <- getStates(this, t);
      } else {
        t <- as.integer(t);
      }
      T[,kk] <- t;
    } # for (kk ...)
  }
  if (byNames) {
    T <- getStateNames(this, T);
  }

  T;
})

setMethodS3("hasState", "MatrixRocData", function(this, state, ...) {
  byNames <- is.character(state);

  T <- getTruth(this, byNames=byNames);
  res <- logical(length(T));
  dim(res) <- dim(T);

  keep <- is.na(T);
  res[keep] <- FALSE;
  keep <- !keep;
  res[keep] <- (T[keep] == state);

  res;
})


setMethodS3("getStatesMap", "MatrixRocData", function(this, ...) {
  statesMap <- this$.statesMap;
  if (is.null(statesMap)) {
    statesMap <- c(NEUTRAL=0, LOSS=-1, GAIN=+1);
  }
  statesMap;
})


setMethodS3("getStateNames", "MatrixRocData", function(this, states) {
  statesMap <- getStatesMap(this);
  idxs <- match(states, statesMap);
  names <- names(statesMap)[idxs];
  dim(names) <- dim(states);
  names;
}, protected=TRUE)

setMethodS3("getStates", "MatrixRocData", function(this, names) {
  statesMap <- getStatesMap(this);
  idxs <- match(names, names(statesMap));
  states <- statesMap[idxs];
  dim(states) <- dim(names);
  states;
}, protected=TRUE)



setMethodS3("plotTracks", "MatrixRocData", function(this, ..., xlim=NULL, ylim=c(-5,5), pch=19, xlab="Position", ylab="Signal") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  colorFcn <- function(states, ...) {
    col <- rep("gray", length(states));
    col[states == 0] <- "black";
    col[states  < 0] <- "red";
    col[states  > 0] <- "blue";
    col;
  } # colorFcn()


  Y <- getData(this);
  T <- getTruth(this, byNames=FALSE);
  x <- seq(length=nrow(Y));

  if (is.null(xlim)) {
    xlim <- range(x, na.rm=TRUE);
  }

  iis <- seq(length=ncol(Y));
  layout(matrix(iis, ncol=1));
  par(mar=c(3,3,1,1)+0.1, pch=pch);
  for (ii in iis) {
    plot(NA, xlim=xlim, ylim=ylim);
    points(x, Y[,ii], col=colorFcn(T[,ii]));
  }
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   BEGIN: Methods that requires matrix data
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("extractSubset", "MatrixRocData", function(this, rows=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- dim(this);
  # Sanity check
  if (is.null(dim)) {
    throw("Cannot extract data by rows. The data has no dimension.");
  }

  # Argument 'rows':
  if (!is.null(rows)) {
    rows <- Arguments$getIndices(rows, range=c(1, dim[1]));
  }

  # Nothing to do?
  if (is.null(rows))
    return(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return a cloned subset
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- clone(this);
  clearCache(res);
  fields <- getInternalRocFields(this);
  for (field in fields) {
    values <- res[[field]];
    if (!identical(length(values), dim[2])) {
      values <- values[rows,,drop=FALSE];
    }
    res[[field]] <- values;
  }

  res;
})




setMethodS3("extractSmoothRocData", "MatrixRocData", function(this, h, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(2, Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Sanity check
  dim <- dim(this);
  if (is.null(dim)) {
    throw("Cannot extract data by rows. The data has no dimension.");
  }

  verbose && enter(verbose, "Creating smoothed data set");
  verbose && cat(verbose, "Amount of smoothing: ", h);

## CACHING is slower than generating from scratch. /HB 2008-07-25
##  verbose && enter(verbose, "Checking for cached results");
##  key <- list(method="extractSmoothRocData", class=class(this)[1], checksum=getChecksum(this), h=h);
##  dirs <- c("MatrixRocData");
##  res <- loadCache(key=key, dirs=dirs);
##  if (!force && !is.null(res)) {
##    verbose && cat(verbose, "Found cached results!");
##    verbose && exit(verbose);
##    verbose && exit(verbose);
##    return(res);
##  }
##  verbose && exit(verbose);

  verbose && enter(verbose, "Cloning existing ", class(this)[1], " object");
  res <- clone(this);
  clearCache(res);
  verbose && exit(verbose);

  verbose && enter(verbose, "Smoothing applicable fields");
  idxs <- NULL;
  fields <- getInternalRocFields(res);
  for (ff in seq(along=fields)) {
    field <- fields[[ff]];
    verbose && enter(verbose, sprintf("Field %d ('%s') of %d", ff, gsub(".", "", field, fixed=TRUE), length(fields)));
    values <- res[[field]];
    verbose && str(verbose, values);

    # Smooth field?
    if (is.matrix(values)) {
      verbose && enter(verbose, "Smoothing");

      if (is.null(idxs)) {
        verbose && enter(verbose, "Generating block-averaging map");
        idxs <- getBlockAverageMap(n=dim[1], h=h);
        hApprox <- attr(idxs, "hApprox");
        verbose && str(verbose, idxs);
        res$smoothingParameters <- list(h=h, hApprox=hApprox);
        verbose && exit(verbose);
      }

      values <- blockAvg(values, idxs);
      verbose && str(verbose, values);

      res[[field]] <- values;

      verbose && exit(verbose);
    }
    rm(values);
    verbose && exit(verbose);
  } # for (ff ...)
  verbose && exit(verbose);

##  # Cache results?
##  if (cache)
##    saveCache(res, key=key, dirs=dirs);

  verbose && exit(verbose);

  res;
})


setMethodS3("scanTpAtFp", "MatrixRocData", function(this, fpRate, hs=seq(from=1, to=10, by=0.1), ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating the TP rate at a given FP rate for various amounts of smoothing");
  verbose && cat(verbose, "Smooth levels:");
  verbose && print(verbose, hs);

  verbose && enter(verbose, "Loading earlier estimates from file cache");
  key <- list(method="extractSmoothRocData", class=class(this)[1], checksum=getChecksum(this), fpRate=fpRate);
  dirs <- c("MatrixRocData", "scanTpAtFp");
  fit <- loadCache(key=key, dirs=dirs);
  verbose && exit(verbose);

#  # Merge with memory cache  
#  fit <- this$.scanTpAtFpFit;

  if (force)
    fit <- NULL;

  data <- getData(this, raw=TRUE);
  truth <- getTruth(this, raw=TRUE);
  if (!is.matrix(truth))
    truth <- matrix(truth, nrow=nrow(data), ncol=ncol(data), byrow=TRUE);

  fit <- scanRocTpAtFp(truth=truth, data=data, fpRate=fpRate, hs=hs, fit=fit, verbose=verbose);

#  this$.scanTpAtFpFit <- fit;

  # Cache results?
  if (cache)
    saveCache(fit, key=key, dirs=dirs);

  verbose && enter(verbose, "Extract smooth levels of interest");
  keep <- whichVector(fit[,"h"] %in% hs);
  fit <- fit[keep,,drop=FALSE];
  verbose && print(verbose, fit);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
})



setMethodS3("findTpAtFpPerRow", "MatrixRocData", function(this, fpRate=0.01, rows=NULL, ..., skip=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- dim(this);
  if (is.null(dim)) {
    throw("Cannot fit TP rate per data row. Data is not specified as a matrix: ", dim);
  }

  # Argument 'rows':
  if (!is.null(rows)) {
    rows <- Arguments$getIndices(rows, range=c(1, dim[1]));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating the TP rate per data row");
  verbose && cat(verbose, "Rows:");
  verbose && str(verbose, rows);

  # Load cached results
  cache <- loadCachedResult(this, method="findTpAtFpPerRow", fpRate=fpRate, rows=NULL, skip=skip);

  # Get existing estimates
  if (is.null(this$.tpDensityList))
    this$.tpDensityList <- list();
  fpRateStr <- sprintf("fpRate=%f", fpRate);
  tpDensity <- this$.tpDensityList[[fpRateStr]];
  if (is.null(tpDensity)) {
    if (force || is.null(cache$result)) {
      tpDensity <- rep(-1, dim[1]);
    } else {
      tpDensity <- cache$result;
    }
    this$.tpDensityList[[fpRateStr]] <- tpDensity;
  }

  # Keep subset of interest?
  if (!is.null(rows))
    tpDensity <- tpDensity[rows];

  verbose && cat(verbose, "Existing TP-rate estimates:");
  verbose && str(verbose, tpDensity);
  verbose && summary(verbose, tpDensity);

  # Identify data rows to update
  if (force) {
    rrs <- seq(along=tpDensity);
  } else {
    rrs <- whichVector(is.finite(tpDensity) & tpDensity < 0);
  }

  # Nothing more to do?
  n <- length(rrs);
  if (n == 0) {
    verbose && cat(verbose, "Already done!");
    verbose && exit(verbose);
    return(tpDensity);
  }

  # Extract truth and data
  verbose && enter(verbose, "Extracting truth and data");
  truth <- getTruth(this, ordered=FALSE, complete=FALSE, rows=rows);
  dimnames(truth) <- NULL;
  verbose && cat(verbose, "Truth:");
  verbose && str(verbose, truth);
  data <- getData(this, ordered=FALSE, complete=FALSE, rows=rows);
  dimnames(data) <- NULL;
  verbose && cat(verbose, "Data:");
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # For each data row
  verbose && enter(verbose, "Processing all data rows");
  verbose && cat(verbose, "Number of data rows to do: ", n);
  for (kk in seq(along=rrs)) {
    rr <- rrs[kk];

    # Store every 100th iteration
    if (kk %% 500 == 0) {
      if (is.null(rows)) {
        this$.tpDensityList[[fpRateStr]] <- tpDensity;
      } else {
        this$.tpDensityList[[fpRateStr]][rows] <- tpDensity;
      }
      verbose && printf(verbose, "Data row #%d of %d\n", kk, n);
    }

    t <- truth[rr,];
    d <- data[rr,];
    fit <- findRocTpAtFp(truth=t, data=d, fpRate=fpRate, hasNAs=TRUE, isOrdered=FALSE, .checkArgs=FALSE);
#    print(data.frame(truth=t, data=d));
#    print(fit);
    tpRate <- fit[["tpRateEst"]]; 
    tpDensity[rr] <- tpRate;
  } # for (kk in ...)
  verbose && exit(verbose);

  verbose && enter(verbose, "Storing estimates in memory cache");
  # Store estimates
  if (is.null(rows)) {
    this$.tpDensityList[[fpRateStr]] <- tpDensity;
  } else {
    this$.tpDensityList[[fpRateStr]][rows] <- tpDensity;
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return updates
  res <- tpDensity;

##  verbose && str(verbose, tpDensity);
##  verbose && summary(verbose, tpDensity);

  # Cache result?
  doCache <- (is.null(rows) || length(rows) > 1000);
  doCache <- TRUE;
  if (doCache) {
    if (!is.null(rows)) {
      tpDensity <- rep(-1, dim[1]);
      tpDensity[rows] <- res;
    }

##    verbose && str(verbose, tpDensity);
##    verbose && summary(verbose, tpDensity);

    if (is.null(cache$key)) {
      cache <- loadCachedResult(this, method="findTpAtFpPerRow", fpRate=fpRate, rows=NULL, skip=FALSE);
    }

##    verbose && str(verbose, tpDensity);
##    verbose && summary(verbose, tpDensity);

    saveCache(tpDensity, key=cache$key, dirs=cache$dirs);
  }

  res;
}) # findTpAtFpPerRow()



setMethodS3("drawTpRateDensity", "MatrixRocData", function(this, fpRate=0.01, rows=NULL, ..., adjust=1, lwd=2, col=par("col"), cex=0.7, annotate=!add, add=FALSE, force=FALSE, verbose=FALSE) {
  require("aroma.core") || throw("Package not loaded: aroma.core"); # stext()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tpRates <- findTpAtFpPerRow(this, fpRate=fpRate, rows=rows, force=force, verbose=verbose);
  n <- length(tpRates);

  # Filter out non-finite estimates
  keep <- is.finite(tpRates);
  tpRates <- tpRates[keep];
  rm(keep);

  callRate <- length(tpRates)/n;
  d <- density(tpRates, from=0, to=1, adjust=adjust);

  if (add) {
    lines(d, lwd=lwd, col=col, ...);
  } else {
    tpLab <- "TP rate";
    ylab <- "Density (integrates to one)";
    plot(d, lwd=lwd, col=col, main="", xlab=tpLab, ylab=ylab, ...);
  }

  if (annotate) {
    text <- sprintf("Call rate: %.3f%% (n=%d)", 100*callRate, n);
    stext(side=3, pos=1, text, cex=cex);
  }

  invisible(d);
})



setMethodS3("plotTpRateDensity", "MatrixRocData", function(this, fpRate=0.01, rows=NULL, ..., col=par("col"), cex=0.7, breaks=NULL, offset=0, width=1, xlab="TP rate", ylab="Count", annotate=!add, add=FALSE, force=FALSE, verbose=FALSE) {
  require("R.basic") || throw("Package not loaded: R.basic"); # hist()
  require("aroma.core") || throw("Package not loaded: aroma.core"); # stext()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'breaks':
  if (!is.null(breaks)) {
    breaks <- Arguments$getDoubles(breaks);
  }

  # Argument 'offset':
  offset <- Arguments$getDouble(offset, range=c(0,Inf));

  # Argument 'width':
  width <- Arguments$getDouble(width, range=c(0,Inf));

  tpRates <- findTpAtFpPerRow(this, fpRate=fpRate, rows=rows, 
                                             force=force, verbose=verbose);
  n <- length(tpRates);

  # Filter out non-finite estimates
  keep <- is.finite(tpRates);
  tpRates <- tpRates[keep];
  rm(keep);

  res <- hist(tpRates, breaks=breaks, offset=offset, width=width, 
              freq=TRUE, col=col, border=NA, add=add, 
              xlab=xlab, ylab=ylab, main="");

  if (annotate) {
    callRate <- length(tpRates)/n;
    text <- sprintf("Call rate: %.3f%% (n=%d)", 100*callRate, n);
    stext(side=3, pos=1, text, cex=cex);
  }

  invisible(res);
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   END: Methods that requires matrix data
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



############################################################################
# HISTORY:
# 2009-02-01
# o Imported to aroma.cn.eval.
# 2008-07-25
# o Now clearCache() clears all cached results.
# o Added extractSmoothedRocData().  This assumes that the rows in the 
#   data matrix is ordered along the genome.
# 2008-07-24
# o Added fast findTpAtFpPerRow().
# o Added fast findTpAtFp().
# 2008-07-23
# o Created RocData class.
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
