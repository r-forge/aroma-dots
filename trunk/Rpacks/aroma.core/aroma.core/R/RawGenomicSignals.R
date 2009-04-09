###########################################################################/**
# @RdocClass RawGenomicSignals
#
# @title "The RawGenomicSignals class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of length J specifying the signal
#     at each loci.}
#   \item{x}{A (optional) @numeric @vector of length J specifying the 
#     position of each loci.}
#   \item{chromosome}{An (optional) @integer specifying the chromosome for
#     these genomic signals.}
#   \item{name}{An (optional) @characte string specifying the sample name.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("RawGenomicSignals", function(y=NULL, x=NULL, chromosome=NA, name=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  if (!is.null(y)) {
    if (inherits(y, "RawGenomicSignals")) {
      object <- y;
      y <- object$y;
      x <- object$x;
      chromosome <- object$chromosome;
      name <- object$name;
      rm(object);
    }

    if (!is.vector(y)) {
      throw("Argument 'y' must be a vector: ", mode(y)[1]);
    }

    if (!is.numeric(y)) {
      throw("Argument 'y' must be a numeric: ", class(y)[1]);
    }
  }

  if (!is.null(x)) {
    if (!is.vector(x)) {
      throw("Argument 'x' must be a vector: ", mode(x)[1]);
    }

    if (!is.numeric(x)) {
      throw("Argument 'x' must be a numeric: ", class(x)[1]);
    }

    n <- length(y);
    if (length(x) != n) {
      throw("Argument 'x' and 'y' are of different lengths: ", length(x), " != ", n);
    }
  }

  # Argument 'chromosome':
  if (!is.na(chromosome)) {
    chromosome <- Arguments$getIndex(chromosome);
  }


  extend(Object(), "RawGenomicSignals", 
    y = y,
    x = x,
    chromosome = chromosome,
    .name = name
  )
})


setMethodS3("as.character", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  name <- getName(this);
  if (is.null(name)) name <- "";
  s <- c(s, sprintf("Name: %s", as.character(name)));
  s <- c(s, sprintf("Chromosome: %d", getChromosome(this)));
  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(this)));
  fields <- getLociFields(this);
  fields <- sapply(fields, FUN=function(field) {
    values <- this[[field]];
    mode <- mode(values);
    sprintf("%s [%dx%s]", field, length(values), mode);
  })
  fields <- paste(fields, collapse=", ");
  s <- c(s, sprintf("Loci fields: %s", fields));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE) 


setMethodS3("nbrOfLoci", "RawGenomicSignals", function(this, na.rm=FALSE, ...) {
  y <- getSignals(this);
  if (na.rm) {
    y <- y[is.finite(y)];
  }
  length(y);
})

setMethodS3("getPositions", "RawGenomicSignals", function(this, ...) {
  x <- this$x;
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci(this));
  }
  x;
})


setMethodS3("getChromosome", "RawGenomicSignals", function(this, ...) {
  chr <- this$chromosome;
  if (is.null(chr))
    chr <- NA;
  chr <- as.integer(chr);
  chr;
})


setMethodS3("getSignals", "RawGenomicSignals", function(this, ...) {
  this$y;
})


setMethodS3("getName", "RawGenomicSignals", function(this, ...) {
  this$.name;
})


setMethodS3("as.data.frame", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  data.frame(x=this$x, y=this$y);
})

setMethodS3("summary", "RawGenomicSignals", function(object, ...) {
  # To please R CMD check
  this <- object;

  summary(as.data.frame(this));
})


setMethodS3("getLociFields", "RawGenomicSignals", function(this, ...) {
  c("x", "y");
})

setMethodS3("sort", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  res <- clone(this);
  x <- getPositions(res);
  o <- order(x);
  rm(x);
  for (field in getLociFields(res)) {
    res[[field]] <- res[[field]][o];
  }
  res;
})


setMethodS3("getXY", "RawGenomicSignals", function(this, sort=TRUE, ...) {
  xy <- data.frame(x=this$x, y=this$y);
  if (sort)
    xy <- xy[order(xy$x),,drop=FALSE];
  xy;
})



setMethodS3("extractSubset", "RawGenomicSignals", function(this, subset, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':
  subset <- Arguments$getIndices(subset, range=c(1, nbrOfLoci(this)));

  y <- getSignals(this);
  x <- getPositions(this);
  y <- y[subset];
  x <- x[subset];

  res <- clone(this);
  clearCache(res);
  res$x <- x;
  res$y <- y;
  rm(y, x);

  res;
})



setMethodS3("kernelSmoothing", "RawGenomicSignals", function(this, xOut=NULL, ..., verbose=FALSE) {
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
  y <- getSignals(this);
  x <- getPositions(this);

  if (is.null(xOut)) {
    xOut <- x;
  }

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  verbose && enter(verbose, "Kernel smoothing");
  verbose && cat(verbose, "Arguments:");
  args <- list(y=y, x=x, xOut=xOut, ...);
  verbose && str(verbose, args);
  ys <- kernelSmoothing(y=y, x=x, xOut=xOut, ...);
  verbose && str(verbose, ys);
  verbose && exit(verbose);


  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$y <- ys;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # kernelSmoothing()


setMethodS3("gaussianSmoothing", "RawGenomicSignals", function(this, sd=10e3, ...) {
  kernelSmoothing(this, kernel="gaussian", h=sd, ...);
})



setMethodS3("binnedSmoothing", "RawGenomicSignals", function(this, ..., byCount=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'byCount':
  byCount <- Arguments$getLogical(byCount);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  y <- getSignals(this);
  x <- getPositions(this);
  verbose && printf(verbose, "Range of positions: [%d,%d]\n", 
                       min(x,na.rm=TRUE), max(x,na.rm=TRUE));

  if (byCount) {
    verbose && enter(verbose, "Binned smoothing (by count)");
    # Smoothing both y and x.
    Y <- cbind(y=y, x=x);
    xRank <- seq(length=length(x));
    verbose && cat(verbose, "Positions (ranks):");
    verbose && str(verbose, xRank);
    verbose && cat(verbose, "Arguments:");
    args <- list(Y=Y, x=x, ...);
    verbose && str(verbose, args);
    Ys <- colBinnedSmoothing(Y=Y, x=xRank, ..., verbose=less(verbose, 10));
    verbose && str(verbose, Ys);
    xOut <- attr(Ys, "xOut");
    verbose && str(verbose, xOut);
    # The smoothed y:s
    ys <- Ys[,1,drop=TRUE];
    # The smoothed x:s, which becomes the new target positions
    xOut <- Ys[,2,drop=TRUE];
    rm(xRank, Y, Ys);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Binned smoothing (by position)");
    verbose && cat(verbose, "Arguments:");
    args <- list(y=y, x=x, ...);
    verbose && str(verbose, args);
    ys <- binnedSmoothing(y=y, x=x, ..., verbose=less(verbose, 10));
    verbose && str(verbose, ys);
    xOut <- attr(ys, "xOut");
    verbose && exit(verbose);
  }


  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$y <- ys;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # binnedSmoothing()




###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod estimateStandardDeviation
#
# @title "Estimates the standard deviation of the raw Ys"
#
# \description{
#  @get "title" robustly or non-robustly using either a "direct" estimator
#  or a first-order difference estimator.
# }
#
# @synopsis
#
# \arguments{
#   \item{method}{If \code{"diff"}, the estimate is based on the first-order
#     contigous differences of raw Ys. If \code{"direct"}, it is based 
#     directly on the raw Ys.}
#   \item{estimator}{If \code{"mad"}, the robust @see "stats::mad" estimator
#     is used.  If \code{"sd"}, the @see "stats::sd" estimator is used.}
#   \item{na.rm}{If @TRUE, missing values are excluded first.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a non-negative @numeric value.
# }
#
# @author
#
# \seealso{
#   @see "base::diff", @see "stats::sd", and @see "stats::mad".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/########################################################################### 
setMethodS3("estimateStandardDeviation", "RawGenomicSignals", function(this, method=c("diff", "direct"), estimator=c("mad", "sd"), na.rm=TRUE, ...) {
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'estimator':
  estimator <- match.arg(estimator);


  # Get the estimator function
  estimatorFcn <- get(estimator, mode="function");

  y <- getSignals(this);
  if (method == "diff") {
    y <- diff(y);
    sigma <- estimatorFcn(y, na.rm=na.rm)/sqrt(2);
  } else if (method == "direct") {
    sigma <- estimatorFcn(y, na.rm=na.rm);
  }

  sigma;
})



setMethodS3("plot", "RawGenomicSignals", function(x, xlab="Position", ylab="Signal", ylim=c(-3,3), pch=20, xScale=1, yScale=1, ...) {
  # To please R CMD check
  this <- x;

  x <- getPositions(this);
  y <- getSignals(this);
  plot(xScale*x, yScale*y, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, ...);
})


setMethodS3("points", "RawGenomicSignals", function(x, pch=20, xScale=1, yScale=1, ...) {
  # To please R CMD check
  this <- x;

  x <- getPositions(this);
  y <- getSignals(this);
  points(xScale*x, yScale*y, pch=pch, ...);
})

setMethodS3("lines", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  x <- getPositions(this);
  y <- getSignals(this);

  o <- order(x);
  x <- x[o];
  y <- y[o];
  lines(x, y, ...);
})


setMethodS3("xSeq", "RawGenomicSignals", function(this, from=1, to=xMax(this), by=100e3, ...) {
  seq(from=from, to=to, by=by);
})

setMethodS3("xRange", "RawGenomicSignals", function(this, na.rm=TRUE, ...) {
  x <- getPositions(this);
  range(x, na.rm=na.rm);
})

setMethodS3("xMin", "RawGenomicSignals", function(this, ...) {
  xRange(this, ...)[1];
})

setMethodS3("xMax", "RawGenomicSignals", function(this, ...) {
  xRange(this, ...)[2];
})

setMethodS3("signalRange", "RawGenomicSignals", function(this, na.rm=TRUE, ...) {
  y <- getSignals(this);
  range(y, na.rm=na.rm);
})



setMethodS3("extractRawGenomicSignals", "default", abstract=TRUE);



############################################################################
# HISTORY:
# 2009-04-06
# o BUG FIX: binnedSmoothing(..., byCount=TRUE) of RawGenomicSignals would
#   give error "[...] object "ys" not found".
# 2009-02-19
# o Renamed from RawCopyNumbers RawGenomicSignals.
# o Added argument 'byCount' to binnedSmoothing() of RawGenomicSignals.
# 2009-02-17
# o Now RawGenomicSignals() also takes another RawCopyNumbers object as
#   input.
# 2009-02-16
# o Added optional constructor argument 'name'.
# 2009-02-07
# o Added Rdoc comments and example.
# 2008-05-21
# o Added field 'chromosome' (single value).
# 2008-05-17
# o Added abstract default extractCopyNumberRegions().
# o Moved to aroma.core. 
# 2008-03-31
# o Put recently added sd() and mad() into estimateStandardDeviation().
# 2008-03-10
# o Added standard deviation estimator sd() and mad() which my default
#   uses a first-order difference variance estimator.
# 2007-08-22
# o Created.  Need a generic container for holding copy number data and
#   to plot them nicely.
############################################################################
