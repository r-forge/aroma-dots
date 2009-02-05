setConstructorS3("MatrixRocData", function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  extend(RocData(...), "MatrixRocData");
}) # MatrixRocData()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   BEGIN: Methods that requires matrix data
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("extractSubset", "MatrixRocData", function(this, rows=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows':
  dim <- dim(this);
  if (!is.null(rows))
    rows <- Arguments$getIndices(rows, range=c(1, dim[1]));

  # Sanity check
  if (!is.null(rows)) {
    if (is.null(dim))
      throw("Cannot extract data by rows. The data has no dimension.");
  }

  # Nothing to do?
  if (is.null(rows))
    return(this);

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
  if (is.null(dim))
    throw("Cannot extract data by rows. The data has no dimension.");

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
