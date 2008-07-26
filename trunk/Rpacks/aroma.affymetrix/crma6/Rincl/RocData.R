# \arguments{
#   \item{truth}{A @numerical @vector of length N or a @numerical NxK @matrix.}
#   \item{data}{A @numerical object of same dimension as \code{truth}.}
# }
setConstructorS3("RocData", function(truth=NULL, data=NULL, recall=NULL, ...) {
  # Argument 'truth':
  if (inherits(truth, "RocData")) {
    object <- truth;
    truth <- object$.truth;
    data <- object$.data;
  }

  # Argument 'data':
  if (!is.null(data)) {
    if (is.numeric(truth)) {
      dim <- dim(data);

      # Is truth specified by columns (in a matrix)?
      truthByColumns <- (!is.null(dim) && length(truth) == dim[2]);
      
      if (!truthByColumns) {
        # Verify same lengths
        if (length(data) != length(truth)) {
          throw("Argument 'data' and 'truth' are of different lengths: ", 
                                   length(data), " != ", length(truth));
        }

        # Verify same dimensions
        if (!identical(dim(truth), dim)) {
          throw("Argument 'truth' and 'data' have different dimensions: ", 
                                  paste(dim(truth), collapse="x"), " != ", 
                                                paste(dim, collapse="x"));
        }
      }

      # Turn 'truth' into (0,1) variable by re-calling?
      if (!is.null(recall))
        truth <- (truth == recall);

      # Turn into (0,1) integers
      truth <- as.integer(truth);

      # Preserve dimensions, if any.
      if (!truthByColumns)
        dim(truth) <- dim;
    }
  }


  this <- extend(Object(), "RocData",
    .truth=truth,
    .data=data
  );

  this;
}) # RocData()


setMethodS3("clearCache", "RocData", function(this, ...) {
  fields <- c(".checksum", ".order", ".hasNAs", ".isMissing", 
              ".rate", ".tpDensityList");
  for (ff in fields) {
    this[[ff]] <- NULL;
  }
  NextMethod("clearCache", this, ...);
})


setMethodS3("getChecksum", "RocData", function(this, force=FALSE, ..., skip=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  checksum <- this$.checksum;

  if (!skip && (force || is.null(checksum))) {
    fields <- getRocFields(this);
    for (ff in seq(along=fields)) {
      field <- fields[[ff]];
      verbose && enter(verbose, sprintf("Calculating checksum for field %d ('%s') of %d", ff, field, length(fields)));
  
      verbose && enter(verbose, "Getting data");
      values <- getField(this, field=field, raw=TRUE);
      verbose && exit(verbose);

      verbose && enter(verbose, "Generating checksum");
      checksumFF <- digest(list(values));
      rm(values);
      verbose && exit(verbose);
  
      # Cumulative checksum
      checksum <- digest(list(checksum, checksumFF));
  
      verbose && exit(verbose);
    } # for (ff ...)

    this$.checksum <- checksum;
  }

  checksum;
})


setMethodS3("as.character", "RocData", function(this, ...) {
  s <- sprintf("%s: ", class(this)[1]);
  s <- c(s, paste("ROC fields:", paste(getRocFields(this), collapse=", ")));
  s <- c(s, sprintf("Data dimension: %s", paste(dim(this), collapse="x")));
  s <- c(s, sprintf("Truth by columns: %s", isTruthByColumns(this)));
  s <- c(s, sprintf("Number of data points: %d", nbrOfDataPoints(this, complete=FALSE)));
  s <- c(s, sprintf("Call rate: %f%%", 100*getCallRate(this)));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("isTruthByColumns", "RocData", function(this, ...) {
  dim <- dim(this);
  if (is.null(dim))
    return(FALSE);

  field <- getInternalRocFields(this, "truth");
  (length(this[[field]]) == dim[2]);
})

setMethodS3("dim", "RocData", function(this, ...) {
  dim(getData(this, ordered=FALSE, complete=FALSE, ...));
})

setMethodS3("getRocFields", "RocData", function(this, ...) {
  c("truth", "data");
})

setMethodS3("getInternalRocFields", "RocData", function(this, fields=NULL, ...) {
  allFields <- getRocFields(this, ...);
  internalFields <- paste(".", allFields, sep="");
  names(internalFields) <- allFields;

  # Keep certain fields?
  if (!is.null(fields)) {
    internalFields <- internalFields[fields];
  }

  internalFields;
})


setMethodS3("getOrder", "RocData", function(this, force=FALSE, ...) {
  o <- this$.order;
  if (force || is.null(o)) {
    values <- getData(this, ordered=FALSE, complete=FALSE);
    o <- order(values);
    this$.order <- o;
  }
  o;
})

# Returns true if any of the ROC data fields have NAs
setMethodS3("anyMissing", "RocData", function(this, ...) {
  hasNAs <- this$.hasNAs;
  if (hasNAs) {
    for (ff in getInternalRocFields(this)) {
      values <- this[[ff]];
      if (anyMissing(values)) {
        hasNAs <- TRUE;
        break;
      }
    }
    this$.hasNAs <- hasNAs;
  }
  hasNAs;
});


# Returns a @logical @vector indicating if whether the correspond
# ROC tuple has a missing value (in either of the ROC data fields).
# If \code{ordered==TRUE}, the ROC tuples are ordered according 
# to the data field, otherwise not.
setMethodS3("isMissing", "RocData", function(this, ordered=TRUE, force=FALSE, ...) {
  isMissing <- this$.isMissing;
  if (force || is.null(isMissing)) {
    isMissing <- rep(FALSE, times=nbrOfDataPoints(this, complete=FALSE));
    for (ff in getInternalRocFields(this)) {
      values <- this[[ff]];
      isMissing <- isMissing | is.na(values);
    }
#    dim(isMissing) <- dim(this);
    this$.isMissing <- isMissing;
  }

  # Ordered by data values?
  if (ordered)
    isMissing <- isMissing[getOrder(this)];

  isMissing;
})

setMethodS3("nbrOfDataPoints", "RocData", function(this, complete=TRUE, ...) {
  data <- getData(this, ordered=FALSE, complete=complete, ...);
  length(data);
})


setMethodS3("getCallRate", "RocData", function(this, ...) {
  rate <- this$.rate;
  if (is.null(rate)) {
    rate <- nbrOfDataPoints(this, complete=TRUE) / nbrOfDataPoints(this, complete=FALSE);
    this$.rate <- rate;
  }
  rate;
})

setMethodS3("as.data.frame", "RocData", function(x, ...) {
  extractDataFrame(x, ...);
})

setMethodS3("extractDataFrame", "RocData", function(this, fields=NULL, ...) {
  knownFields <- getRocFields(this);
  if (is.null(fields)) {
    fields <- knownFields; 
  } else if (any(!fields %in% knownFields)) {
    missing <- fields[!(fields %in% knownFields)];
    throw("Unkown field(s): ", paste(missing, collapse=", "));
  }

  internalFields <- getInternalRocFields(this, fields=fields);
  # Pre-allocate data frame
  modes <- sapply(internalFields, FUN=function(ff) {
    storage.mode(this[[ff]]);
  });
  colClasses <- modes;
  names(colClasses) <- fields;

  df <- NULL;

  # Extract 
  for (ff in seq(along=fields)) {
    field <- fields[ff];
    values <- getField(this, field=field, ...);
    storage.mode(values) <- modes[ff];
    if (is.null(df)) {
      df <- dataFrame(colClasses=colClasses, nrow=length(values));
    }
    df[[ff]] <- values;
    rm(values);
  }

  df;
})


setMethodS3("getField", "RocData", function(this, field, raw=FALSE, ordered=TRUE, complete=TRUE, rows=NULL, ...) {
  internalField <- getInternalRocFields(this, fields=field);
  values <- this[[internalField]];

  if (raw) {
    return(values);
  }

  # Expand 'truth', if specified as one value per column
  if (field == "truth") {
    dim <- dim(this);
    if (!is.null(dim)) {
      colTruth <- values;
      values <- matrix(colTruth[1], nrow=dim[1], ncol=dim[2]);
      for (cc in seq(along=colTruth)) {
        values[,cc] <- colTruth[cc];
      }
    }
  }

  if (!is.null(rows)) {
    dim <- dim(this);
    if (is.null(dim))
      throw("Cannot extract data by rows. The data has no dimension.");
    values <- values[rows,,drop=FALSE];
  }

  # Order it?
  if (ordered)
    values <- values[getOrder(this)];

  # Exclude non-complete ROC tuples?
  if (complete)
    values <- values[!isMissing(this, ordered=ordered)];

  values;
}, protected=TRUE)


setMethodS3("getTruth", "RocData", function(this, ...) {
  getField(this, field="truth", ...);
})

setMethodS3("getData", "RocData", function(this, ...) {
  getField(this, field="data", ...);
})


setMethodS3("loadCachedResult", "RocData", function(this, method, ..., skip=FALSE, dirs=NULL) {
  # Get checksum for this object
  checksum <- getChecksum(this, skip=skip);
  if (is.null(checksum))
    return(NULL);

  key <- list(method=method, class=class(this), checksum=checksum, ...);

  if (is.null(dirs))
    dirs <- c(class(this)[1], method);
  
  res <- loadCache(key=key, dirs=dirs);
  list(result=res, key=key, dirs=dirs);
})


setMethodS3("fit", "RocData", function(this, ..., skip=FALSE, force=FALSE) {
  # Check file cache for results
  cache <- loadCachedResult(this, method="fit", complete=TRUE, skip=skip);
  if (!force && !is.null(cache$result))
    return(cache$result);

  truth <- getTruth(this, ordered=TRUE, complete=TRUE);
  data <- getData(this, ordered=TRUE, complete=TRUE);
#str(list(truth=truth, data=data));
  fit <- fitRocLite(truth=truth, data=data, hasNAs=FALSE, ncuts=1200, isOrdered=TRUE);
  res <- RocCurve(roc=fit$roc, cuts=fit$cuts);

  # Cache result
  saveCache(res, key=cache$key, dirs=cache$dirs);

  res;
})


setMethodS3("findTpAtFp", "RocData", function(this, ..., skip=FALSE) {
  # Check file cache for results
  cache <- loadCachedResult(this, method="findTpAtFp", ..., skip=skip);
  if (!force && !is.null(cache$result))
    return(cache$result);    

  truth <- getTruth(this, ordered=TRUE, complete=TRUE);
  data <- getData(this, ordered=TRUE, complete=TRUE);
  res <- findTpAtFpLite(truth=truth, data=data, hasNAs=FALSE, isOrdered=TRUE, ...);

  # Cache result
  saveCache(res, key=cache$key, dirs=cache$dirs);

  res;
}) # findTpAtFp()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   BEGIN: Methods that requires matrix data
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("extractSubset", "RocData", function(this, rows=NULL, ...) {
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




setMethodS3("extractSmoothRocData", "RocData", function(this, h, ..., verbose=FALSE) {
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
##  dirs <- c("RocData");
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



setMethodS3("findTpAtFpPerRow", "RocData", function(this, fpRate=0.01, rows=NULL, ..., skip=FALSE, force=FALSE, verbose=FALSE) {
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
    if (kk %% 100 == 0) {
      if (is.null(rows)) {
        this$.tpDensityList[[fpRateStr]] <- tpDensity;
      } else {
        this$.tpDensityList[[fpRateStr]][rows] <- tpDensity;
      }
      verbose && printf(verbose, "Data row #%d of %d\n", kk, n);
    }

    t <- truth[rr,];
    d <- data[rr,];
    fit <- findTpAtFpLite(truth=t, data=d, fpRate=fpRate, hasNAs=TRUE, isOrdered=FALSE, .checkArgs=FALSE);
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



setMethodS3("plotTpRateDensity", "RocData", function(this, fpRate=0.01, rows=NULL, ..., lwd=2, col=par("col"), cex=0.7, annotate=!add, add=FALSE, force=FALSE, verbose=FALSE) {
  tpRates <- findTpAtFpPerRow(this, fpRate=fpRate, rows=rows, force=force, verbose=verbose);
  n <- length(tpRates);

  # Filter out non-finite estimates
  keep <- is.finite(tpRates);
  tpRates <- tpRates[keep];
  rm(keep);

  callRate <- length(tpRates)/n;
  d <- density(tpRates, from=0, to=1);

  if (add) {
    lines(d, lwd=lwd, col=col);
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   END: Methods that requires matrix data
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



############################################################################
# HISTORY:
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
