###########################################################################/**
# @RdocClass SegmentedCopyNumbers
#
# @title "The SegmentedCopyNumbers class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RawCopyNumbers".}
#   \item{states}{A @function returning the copy-number states given a
#     @vector of locus positions.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/SegmentedCopyNumbers.Rex"
#
# @author
#*/########################################################################### 
setConstructorS3("SegmentedCopyNumbers", function(..., states=NULL) {
  if (!is.null(states)) {
    if (is.function(states)) {
    } else {
      throw("Argument 'states' must be a function: ", mode(states));
    }
  }

  extend(RawCopyNumbers(...), "SegmentedCopyNumbers", 
    states = states
  )
})

setMethodS3("getStates", "SegmentedCopyNumbers", function(this, x=getPositions(this), ...) {
  states <- this$states;

  if (is.function(states)) {
    fcn <- states;
    chromosome <- getChromosome(this);
    states <- fcn(x, chromosome=chromosome, ...);
    storage.mode(states) <- "integer";
  }

  states;
})



setMethodS3("kernelSmoothingByState", "SegmentedCopyNumbers", function(this, xOut=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'xOut':
  if (!is.null(xOut)) {
    xOut <- Arguments$getDoubles(xOut);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  x <- getPositions(this);

  if (is.null(xOut)) {
    xOut <- x;
  }
  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  naValue <- as.double(NA);
  yOut <- rep(naValue, length(xOut));

  y <- getCNs(this);
  states <- getStates(this);

  uniqueStates <- unique(states);
  uniqueStates <- sort(uniqueStates, na.last=TRUE);
  verbose && cat(verbose, "Unique states:");
  verbose && str(verbose, uniqueStates);

  # Identify states of target loci
  statesOut <- getStates(this, x=xOut);

  for (ss in seq(along=uniqueStates)) {
    state <- uniqueStates[ss];
    verbose && enter(verbose, sprintf("State #%d ('%d') of %d", 
                                      ss, state, length(uniqueStates)));

    # Identifying loci with this state
    if (is.na(state)) {
      keep <- is.na(states);
    } else {
      keep <- (states == state);
    }
    keep <- whichVector(keep);
    statesSS <- states[keep];
    ySS <- y[keep];
    xSS <- x[keep];

    # Identify target loci with this state
    if (is.na(state)) {
      keep <- is.na(statesOut);
    } else {
      keep <- (statesOut == state);
    }
    keep <- whichVector(keep);
    xOutSS <- xOut[keep];

    verbose && enter(verbose, "Kernel smoothing");
    verbose && cat(verbose, "Arguments:");
    args <- list(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, args);
    yOutSS <- kernelSmoothing(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, yOutSS);
    verbose && exit(verbose);
      
    yOut[keep] <- yOutSS;
    verbose && exit(verbose);
  } # for (ss ...)
  verbose && str(verbose, yOut);

  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$cn <- yOut;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # kernelSmoothingByState()



setMethodS3("binnedSmoothingByState", "SegmentedCopyNumbers", function(this, from=xMin(this), to=xMax(this), by=NULL, length.out=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  x <- getPositions(this);
  # Argument 'from' & 'to':
  if (is.null(from)) {
    from <- min(x, na.rm=TRUE);
  } else {
    from <- Arguments$getInteger(from);
  }
  if (is.null(to)) {
    to <- max(x, na.rm=TRUE);
  } else {
    to <- Arguments$getInteger(to, range=c(from, Inf));
  }

  # Arguments 'by' & 'length.out':
  if (is.null(by) & is.null(length.out)) {
    throw("Either argument 'by' or 'length.out' needs to be given.");
  }
  if (!is.null(by)) {
    by <- Arguments$getDouble(by, range=c(0,to-from));
  }
  if (!is.null(length.out)) {
    length.out <- Arguments$getInteger(length.out, range=c(1,Inf));
  }
 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Smoothing data set");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate output loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(by)) {
    xOut <- seq(from=from, to=to, by=by);
  } else {
    xOut <- seq(from=from, to=to, length.out=length.out);
  }
  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);
  naValue <- as.double(NA);
  yOut <- rep(naValue, length(xOut));
  verbose && cat(verbose, "yOut:");
  verbose && str(verbose, yOut);

  statesOut <- getStates(this, x=xOut);
  verbose && cat(verbose, "statesOut:");
  verbose && str(verbose, statesOut);

  uniqueStates <- unique(statesOut);
  uniqueStates <- sort(uniqueStates, na.last=TRUE);
  verbose && cat(verbose, "Unique output states:");
  verbose && str(verbose, uniqueStates);


  y <- getCNs(this);
  states <- getStates(this);

  for (ss in seq(along=uniqueStates)) {
    state <- uniqueStates[ss];
    verbose && enter(verbose, sprintf("State #%d ('%d') of %d", 
                                      ss, state, length(uniqueStates)));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identifying loci with this state
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.na(state)) {
      keep <- is.na(states);
    } else {
      keep <- (states == state);
    }
    keep <- whichVector(keep);
    # Nothing to do?
    if (length(keep) == 0) {
      verbose && exit(verbose);
      next;
    }
    ySS <- y[keep];
    xSS <- x[keep];
    verbose && cat(verbose, "ySS:");
    verbose && str(verbose, ySS);
    verbose && cat(verbose, "xSS:");
    verbose && str(verbose, xSS);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify target loci with this state
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.na(state)) {
      keepOut <- is.na(statesOut);
    } else {
      keepOut <- (statesOut == state);
    }
    keepOut <- whichVector(keepOut);
    verbose && summary(verbose, keepOut);
    statesOutSS <- statesOut[keepOut];
    verbose && cat(verbose, "statesOutSS:");
    verbose && str(verbose, statesOutSS);
    xOutSS <- xOut[keepOut];
    verbose && cat(verbose, "xOutSS:");
    verbose && str(verbose, xOutSS);


    verbose && enter(verbose, "Binned smoothing");
    verbose && cat(verbose, "Arguments:");
    args <- list(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, args);
    yOutSS <- binnedSmoothing(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, yOutSS);
    verbose && exit(verbose);

    yOut[keepOut] <- yOutSS;
    verbose && exit(verbose);
  } # for (ss ...)
  verbose && str(verbose, yOut);

  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$cn <- yOut;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # binnedSmoothingByState()



setMethodS3("getStateColors", "SegmentedCopyNumbers", function(this, ...) {
  states <- getStates(this);

  # Copy neutral states
  col <- rep("#000000", nbrOfLoci(this));

  # Losses
  col[states < 0] <- "blue";

  # Gains
  col[states > 0] <- "red";

  # Unknown
  col[is.na(states)] <- "#999999";

  col;
})


setMethodS3("plot", "SegmentedCopyNumbers", function(x, ..., col=getStateColors(x)) {
  NextMethod("plot", ..., col=col);
})

setMethodS3("points", "SegmentedCopyNumbers", function(x, ..., col=getStateColors(x)) {
  NextMethod("points", ..., col=col);
})



############################################################################
# HISTORY:
# 2009-02-07
# o Created.
############################################################################
