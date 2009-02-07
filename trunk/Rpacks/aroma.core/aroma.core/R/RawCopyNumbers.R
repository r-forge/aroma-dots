###########################################################################/**
# @RdocClass RawCopyNumbers
#
# @title "The RawCopyNumbers class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{cn}{A @numeric @vector of length J specifying the copy number
#     at each loci.}
#   \item{x}{A (optional) @numeric @vector of length J specifying the 
#     position of each loci.}
#   \item{chromosome}{An (optional) @integer specifying the chromosome for
#     these copy numbers.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/RawCopyNumbers.Rex"
#
# @author
#*/########################################################################### 
setConstructorS3("RawCopyNumbers", function(cn=NULL, x=NULL, chromosome=NA, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cn':
  if (!is.null(cn)) {
    if (!is.vector(cn)) {
      throw("Argument 'cn' must be a vector: ", mode(cn)[1]);
    }

    if (!is.numeric(cn)) {
      throw("Argument 'cn' must be a numeric: ", class(cn)[1]);
    }
  }

  if (!is.null(x)) {
    if (!is.vector(x)) {
      throw("Argument 'x' must be a vector: ", mode(x)[1]);
    }

    if (!is.numeric(x)) {
      throw("Argument 'x' must be a numeric: ", class(x)[1]);
    }

    n <- length(cn);
    if (length(x) != n) {
      throw("Argument 'x' and 'cn' are of different lengths: ", length(x), " != ", n);
    }
  }

  # Argument 'chromosome':
  if (!is.na(chromosome)) {
    chromosome <- Arguments$getIndex(chromosome);
  }


  extend(Object(), "RawCopyNumbers", 
    cn = cn,
    chromosome = chromosome,
    x = x
  )
})


setMethodS3("as.character", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
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


setMethodS3("nbrOfLoci", "RawCopyNumbers", function(this, na.rm=FALSE, ...) {
  cn <- this$cn;
  if (na.rm) {
    cn <- cn[is.finite(cn)];
  }
  length(cn);
})

setMethodS3("getPhysicalPositions", "RawCopyNumbers", function(this, ...) {
  getPositions(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("getPositions", "RawCopyNumbers", function(this, ...) {
  x <- this$x;
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci(this));
  }
  x;
})


setMethodS3("getChromosome", "RawCopyNumbers", function(this, ...) {
  chr <- this$chromosome;
  if (is.null(chr))
    chr <- NA;
  chr <- as.integer(chr);
  chr;
})


setMethodS3("getCNs", "RawCopyNumbers", function(this, ...) {
  this$cn;
})

setMethodS3("as.data.frame", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  data.frame(x=this$x, cn=this$cn);
})

setMethodS3("summary", "RawCopyNumbers", function(object, ...) {
  # To please R CMD check
  this <- object;

  summary(as.data.frame(this));
})


setMethodS3("getLociFields", "RawCopyNumbers", function(this, ...) {
  c("cn", "x");
})

setMethodS3("sort", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  res <- clone(this);
  o <- order(res$x);
  for (field in getLociFields(res)) {
    res[[field]] <- res[[field]][o];
  }
  res;
})


setMethodS3("getXY", "RawCopyNumbers", function(this, sort=TRUE, ...) {
  xy <- data.frame(x=this$x, y=this$cn);
  if (sort)
    xy <- xy[order(xy$x),];
  xy;
})



setMethodS3("extractSubset", "RawCopyNumbers", function(this, subset, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':
  subset <- Arguments$getIndices(subset, range=c(1, nbrOfLoci(this)));

  cn <- getCNs(this);
  x <- getPositions(this);
  cn <- cn[subset];
  x <- x[subset];

  res <- clone(this);
  clearCache(res);
  res$cn <- cn;
  res$x <- x;
  rm(cn, x);

  res;
})



setMethodS3("kernelSmoothing", "RawCopyNumbers", function(this, xOut=NULL, ..., verbose=FALSE) {
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
  y <- getCNs(this);
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
  res$cn <- ys;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # kernelSmoothing()


setMethodS3("gaussianSmoothing", "RawCopyNumbers", function(this, sd=10e3, ...) {
  kernelSmoothing(this, kernel="gaussian", h=sd, ...);
})



setMethodS3("binnedSmoothing", "RawCopyNumbers", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  y <- getCNs(this);
  x <- getPositions(this);

  verbose && enter(verbose, "Binned smoothing");
  verbose && cat(verbose, "Arguments:");
  args <- list(y=y, x=x, ...);
  verbose && str(verbose, args);
  ys <- binnedSmoothing(y=y, x=x, ...);
  verbose && str(verbose, ys);
  xOut <- attr(ys, "xOut");
  verbose && exit(verbose);


  verbose && enter(verbose, "Creating result object");
  res <- clone(this);
  clearCache(res);
  res$cn <- ys;
  res$x <- xOut;
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # binnedSmoothing()




###########################################################################/**
# @set "class=RawCopyNumbers"
# @RdocMethod estimateStandardDeviation
#
# @title "Estimates the standard deviation of the raw CNs"
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
#     contigous differences of raw CNs. If \code{"direct"}, it is based 
#     directly on the raw CNs.}
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
setMethodS3("estimateStandardDeviation", "RawCopyNumbers", function(this, method=c("diff", "direct"), estimator=c("mad", "sd"), na.rm=TRUE, ...) {
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'estimator':
  estimator <- match.arg(estimator);


  # Get the estimator function
  estimatorFcn <- get(estimator, mode="function");

  if (method == "diff") {
    value <- diff(this$cn);
    sigma <- estimatorFcn(value, na.rm=na.rm)/sqrt(2);
  } else if (method == "direct") {
    value <- this$cn;
    sigma <- estimatorFcn(value, na.rm=na.rm);
  }

  sigma;
})



setMethodS3("plot", "RawCopyNumbers", function(x, xlab="Position", ylab="Copy number", ylim=c(-3,3), pch=20, xScale=1, yScale=1, ...) {
  # To please R CMD check
  this <- x;

  plot(xScale*this$x, yScale*this$cn, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, ...);
})


setMethodS3("points", "RawCopyNumbers", function(x, pch=20, xScale=1, yScale=1, ...) {
  # To please R CMD check
  this <- x;

  points(xScale*this$x, yScale*this$cn, pch=pch, ...);
})

setMethodS3("lines", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  x <- this$x;
  o <- order(x);
  x <- x[o];
  y <- this$cn[o];
  lines(x, y, ...);
})


setMethodS3("xSeq", "RawCopyNumbers", function(this, from=1, to=xMax(this), by=100e3, ...) {
  seq(from=from, to=to, by=by);
})

setMethodS3("xRange", "RawCopyNumbers", function(this, na.rm=TRUE, ...) {
  values <- this[["x"]];
  range(values, na.rm=na.rm);
})

setMethodS3("xMin", "RawCopyNumbers", function(this, ...) {
  xRange(this, ...)[1];
})

setMethodS3("xMax", "RawCopyNumbers", function(this, ...) {
  xRange(this, ...)[2];
})

setMethodS3("cnRange", "RawCopyNumbers", function(this, na.rm=TRUE, ...) {
  values <- this[["cn"]];
  range(values, na.rm=na.rm);
})



setMethodS3("extractRawCNs", "default", function(...) {
  extractRawCopyNumbers(...);
})

setMethodS3("extractRawCopyNumbers", "default", abstract=TRUE);



############################################################################
# HISTORY:
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
