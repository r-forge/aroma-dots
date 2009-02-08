###########################################################################/**
# @RdocClass CopyNumberRocData
#
# @title "The CopyNumberRocData class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{data}{A @numerical JxK @matrix consisting of L loci in K channels.}
#   \item{truth}{A @function or a @list of K @functions.}
#   \item{positions}{An optional @numeric @vector of J positions.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/CopyNumberRocData.Rex"
#
# @author
#*/########################################################################### 
setConstructorS3("CopyNumberRocData", function(data=NULL, truth=NULL, positions=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'truth':
  if (inherits(data, "CopyNumberRocData")) {
    object <- data;
    data <- object$.data;
    truth <- object$.truth;
    rm(object);
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
    nbrOfLoci <- dim[1];
    nbrOfChannels <- dim[2];

    if (nbrOfChannels == 0) {
      throw("Argument 'data' contains no channels.");
    }

    if (nbrOfLoci < 2) {
      throw("Argument 'data' contains less than two loci: ", nbrOfLoci);
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
    } else {
      throw("Argument 'truth' must be either a function or a list of function: ", mode(truth));
    }

    # Argument 'positions':
    if (!is.null(positions)) {
      positions <- Arguments$getDoubles(positions, length=nbrOfLoci);
    }
  }


  extend(Object(), "CopyNumberRocData",
    .truth=truth,
    .data=data,
    .positions=positions
  )
}) # CopyNumberRocData()


setMethodS3("as.character", "CopyNumberRocData", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s: ", class(this)[1]);
  s <- c(s, sprintf("Data dimension: %s", paste(dim(this), collapse="x")));
  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(this)));
  s <- c(s, sprintf("Number of channels: %d", nbrOfChannels(this)));
#  s <- c(s, sprintf("Call rate: %f%%", 100*getCallRate(this)));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("dim", "CopyNumberRocData", function(x, ...) {
  # To please R CMD check
  this <- x;

  Y <- getData(this);
  dim(Y);
})

setMethodS3("nbrOfLoci", "CopyNumberRocData", function(this, ...) {
  dim(this)[1];
})

setMethodS3("nbrOfChannels", "CopyNumberRocData", function(this, ...) {
  dim(this)[2];
})


setMethodS3("getData", "CopyNumberRocData", function(this, channels=NULL, ...) {
  # Argument 'channels':
  if (!is.null(channels)) {
    nbrOfChannels <- nbrOfChannels(this);
    channels <- Arguments$getIndices(channels, range=c(1, nbrOfChannels));
  }
  
  res <- this$.data;

  if (!is.null(channels)) {
    res <- res[,channels,drop=FALSE];
  }
  
  res;
})

setMethodS3("setData", "CopyNumberRocData", function(this, Y, ...) {
  oldY <- getData(this);

  # Argument 'Y':
  Y <- Arguments$getDoubles(Y);

  if (any(dim(Y) != dim(oldY))) {
    throw("Cannot set new data. Incompatible dimensions: ", paste(dim(Y), collapse="x"), " != ", paste(dim(oldY), collapse="x"));
  }

  this$.data <- Y;

  invisible(this);
})

setMethodS3("getPositions", "CopyNumberRocData", function(this, channels=NULL, ...) {
  # Argument 'channels':
  if (!is.null(channels)) {
    nbrOfChannels <- nbrOfChannels(this);
    channels <- Arguments$getIndices(channels, range=c(1, nbrOfChannels));
  }
  
  res <- this$.positions;

  if (is.null(res)) {
    nbrOfPositions <- nbrOfLoci(this);
    res <- seq(length=nbrOfPositions);
  }

  if (is.matrix(res)) {
    if (!is.null(channels)) {
      res <- res[,channels,drop=FALSE];
    }
  }

  res;
})

setMethodS3("getTruth", "CopyNumberRocData", function(this, channels=NULL, ...) {
  nbrOfLoci <- nbrOfLoci(this);
  nbrOfChannels <- nbrOfChannels(this);

  if (is.null(channels)) {
    channels <- seq(length=nbrOfChannels);
  } else {
    channels <- Arguments$getIndices(channels, range=c(1, nbrOfChannels));
  }

  # Allocate results
  T <- matrix(0L, nrow=nbrOfLoci, ncol=length(channels));
  for (cc in channels) {
    x <- getPositions(this, channels=cc);
    T[,cc] <- getState(this, x=x, channel=cc, ...);
  } # for (cc ...)

  T;
})


setMethodS3("getState", "CopyNumberRocData", function(this, x=NULL, channel, ...) {
  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getDoubles(x);
  }
  
  # Argument 'channel':
  channel <- Arguments$getInteger(channel, range=c(1, nbrOfChannels(this)));

  fcn <- this$.truth;
  if (is.list(fcn)) {
    fcn <- fcn[[channel]];
  }

  if (is.null(x)) {
    x <- getPositions(this, channels=channel);
  }
  t <- fcn(x, channel=channel, ...);
  t <- as.integer(t);

  t;
})


setMethodS3("hasState", "CopyNumberRocData", function(this, state, ...) {
  # Argument 'state':
  state <- Arguments$getInteger(state);

  T <- getTruth(this);
  res <- logical(length(T));
  dim(res) <- dim(T);

  keep <- is.na(T);
  res[keep] <- FALSE;
  keep <- !keep;
  res[keep] <- (T[keep] == state);

  res;
})



setMethodS3("plotTracks", "CopyNumberRocData", function(this, ..., channels=NULL, xlim=NULL, ylim=c(-5,5), pch=19, xlab="Position", ylab="Signal") {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfChannels <- nbrOfChannels(this);
  if (is.null(channels)) {
    channels <- seq(length=nbrOfChannels);
  } else {
    channels <- Arguments$getIndices(channels, range=c(1, nbrOfChannels));
  }


  if (is.null(xlim)) {
    x <- getPositions(this);
    xlim <- range(x, na.rm=TRUE);
    rm(x);
  }

  layout(matrix(seq(along=channels), ncol=1));
  par(mar=c(3,3,1,1)+0.1, pch=pch);
  for (cc in channels) {
    plot(NA, xlim=xlim, ylim=ylim);
    x <- getPositions(this, channel=cc);
    Y <- getData(this, channel=cc);
    T <- getState(this, x=x, channel=cc);
    col <- colorFcn(T);
    points(x, Y, col=col);
  }
})


setMethodS3("extractSubset", "CopyNumberRocData", function(this, rows, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows':
  rows <- Arguments$getIndices(rows, range=c(1, nbrOfLoci(this)));

  data <- getData(this);
  positions <- getPositions(this);
  data <- data[rows,,drop=FALSE];
  if (is.matrix(positions)) {
    positions <- positions[rows,,drop=FALSE];
  } else {
    positions <- positions[rows];
  }

  res <- clone(this);
  clearCache(res);
  res$.data <- data;
  res$.positions <- positions;
  rm(data, positions);

  res;
})



setMethodS3("smooth", "CopyNumberRocData", function(this, xOut=NULL, ..., verbose=FALSE) {
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
  Y <- getData(this);
  x <- getPositions(this);

  if (is.null(xOut)) {
    xOut <- x;
  }

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  verbose && enter(verbose, "Kernel smoothing");
  verbose && cat(verbose, "Arguments:");
  args <- list(Y=Y, x=x, xOut=xOut, ...);
  verbose && str(verbose, args);
  Ys <- colKernelSmoothing(Y=Y, x=x, xOut=xOut, ...);
  verbose && str(verbose, Ys);
  verbose && exit(verbose);


  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$.data <- Ys;
  res$.positions <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # smooth()




setMethodS3("smoothByState", "CopyNumberRocData", function(this, xOut=NULL, ..., verbose=FALSE) {
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

  nbrOfChannels <- nbrOfChannels(this);

  naValue <- as.double(NA);
  Ys <- matrix(naValue, nrow=length(xOut), ncol=nbrOfChannels);

  for (cc in seq(length=nbrOfChannels)) {
    verbose && enter(verbose, sprintf("Channel #%d of %d", cc, nbrOfChannels));

    Ycc <- getData(this, channels=cc);
    Tcc <- getState(this, x=x, channel=cc);
    x <- getPositions(this, channels=cc);

    states <- unique(as.vector(Tcc));
    states <- sort(states, na.last=FALSE);
    verbose && cat(verbose, "CN states:");
    verbose && str(verbose, states);

    # Identify CN states of target loci
    TOut <- getState(this, x=xOut, channel=cc);

    naValue <- as.double(NA);
    Yscc <- rep(naValue, length=length(xOut));

    for (ss in seq(along=states)) {
      state <- states[ss];
      verbose && enter(verbose, sprintf("State #%d ('%d') of %d", 
                                             ss, state, length(states)));

      # Identifying loci with this state
      if (is.na(state)) {
        keep <- is.na(Tcc);
      } else {
        keep <- (Tcc == state);
      }
      keep <- whichVector(keep);
      Tccss <- Tcc[keep];
      Yccss <- Ycc[keep];
      xss <- x[keep];

      # Identify target loci with this state
      if (is.na(state)) {
        keep <- is.na(Tcc);
      } else {
        keep <- (TOut == state);
      }
      keep <- whichVector(keep);
      xOutss <- xOut[keep];

      verbose && enter(verbose, "Kernel smoothing");
      verbose && cat(verbose, "Arguments:");
      args <- list(Y=Yccss, x=xss, xOut=xOutss, ...);
      verbose && str(verbose, args);
      Ysccss <- kernelSmoothing(y=Yccss, x=xss, xOut=xOutss, ...);
      verbose && str(verbose, Ysccss);
      verbose && exit(verbose);
      
      Yscc[keep] <- Ysccss;
      verbose && exit(verbose);
    } # for (ss ...)
    verbose && str(verbose, Yscc);

    Ys[,cc] <- Yscc;

    verbose && exit(verbose);
  } # for (cc ...)

  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$.data <- Ys;
  res$.positions <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # smoothByState()



############################################################################
# HISTORY:
# 2009-02-06
# o Added extractSubset(), smooth(), and smoothByState().
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
