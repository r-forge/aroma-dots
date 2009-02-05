###########################################################################/**
# @RdocClass RocData
#
# @title "The RocData class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{truth}{A @numerical @vector of length N.}
#   \item{data}{A @numerical object of length N.}
#   \item{recall}{(Optional) Unless \code{truth} is given as binary
#     \eqn{\{0,1\}} values, it can be reclassified as such.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Results are memoized, that is, cached to file, in order to speed up
#   subsequent processing.
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("RocData", function(truth=NULL, data=NULL, recall=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


  extend(Object(), "RocData",
    .truth=truth,
    .data=data
  )
}) # RocData()



setMethodS3("clearCache", "RocData", function(this, ...) {
  fields <- c(".checksum", ".order", ".hasNAs", ".isMissing", 
              ".rate", ".tpDensityList", ".scanTpAtFpFit");
  for (ff in fields) {
    this[[ff]] <- NULL;
  }
  NextMethod("clearCache", this, ...);
})


setMethodS3("getChecksum", "RocData", function(this, force=FALSE, ..., skip=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


setMethodS3("as.character", "RocData", function(x, ...) {
  # To please R CMD check
  this <- x;

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
}, private=TRUE)


setMethodS3("dim", "RocData", function(x, ...) {
  # To please R CMD check
  this <- x;

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
      values <- matrix(colTruth, nrow=dim[1], ncol=dim[2], byrow=TRUE);
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
    dirs <- c("aroma.cn.eval", class(this)[1], method);
  
  res <- loadCache(key=key, dirs=dirs);
  list(result=res, key=key, dirs=dirs);
})


setMethodS3("fit", "RocData", function(this, ..., complete=FALSE, skip=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Check file cache for results
  cache <- loadCachedResult(this, method="fit", complete=complete, skip=skip);
  if (!force && !is.null(cache$result))
    return(cache$result);

  verbose && enter(verbose, "Fitting ROC data");

  truth <- getTruth(this, ordered=TRUE, complete=complete);
  data <- getData(this, ordered=TRUE, complete=complete);

  if (!complete) {
    # Treat NAs as FPs
    nas <- whichVector(is.na(data));
    if (length(nas) > 0) {
      verbose && enter(verbose, "Treating NAs as FPs");
      values <- +Inf;
      values <- c(-Inf, Inf);
      values <- rep(values, length.out=length(nas));
      verbose && printf(verbose, "Missing data detected: %d (%.3f%%)\n", length(nas), 100*length(nas)/length(data));
      verbose && cat(verbose, "Summary of imputed missing data:");
      verbose && print(verbose, summary(values > 0));
      data[nas] <- values;
      verbose && exit(verbose);
    }
  }

#str(list(truth=truth, data=data));
  fit <- fitRoc(truth=truth, data=data, hasNAs=FALSE, ncuts=1200, isOrdered=TRUE);

  res <- RocCurve(roc=fit$roc, cuts=fit$cuts);

  # Cache result
  saveCache(res, key=cache$key, dirs=cache$dirs);

  verbose && exit(verbose);

  res;
})


setMethodS3("findTpAtFp", "RocData", function(this, ..., skip=FALSE) {
  # Check file cache for results
  cache <- loadCachedResult(this, method="findTpAtFp", ..., skip=skip);
  if (!force && !is.null(cache$result))
    return(cache$result);    

  truth <- getTruth(this, ordered=TRUE, complete=TRUE);
  data <- getData(this, ordered=TRUE, complete=TRUE);

  res <- findRocTpAtFp(truth=truth, data=data, hasNAs=FALSE, isOrdered=TRUE, ...);

  # Cache result
  saveCache(res, key=cache$key, dirs=cache$dirs);

  res;
}) # findTpAtFp()




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
