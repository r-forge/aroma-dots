###########################################################################/**
# @RdocClass BasePositionNormalization
#
# @title "The BasePositionNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in positioning of
#  A, C, G, and T:s in the probe sequences.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "AbstractProbeSequenceNormalization".}
#   \item{model}{A @character string specifying the model used to fit 
#     the base-count effects.}
#   \item{bootstrap}{If @TRUE, the model fitting is done by bootstrap in
#     order to save memory.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Requirements}{
#   This class requires that an aroma probe sequence file is available
#   for the chip type.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("BasePositionNormalization", function(..., model=c("smooth.spline"), bootstrap=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'bootstrap':
  bootstrap <- Arguments$getLogical(bootstrap);

  if (bootstrap && model != "smooth.spline") {
    throw("Bootstrapping for models other than 'smooth.spline' is not implemented: ", model);
  }


  extend(AbstractProbeSequenceNormalization(...), "BasePositionNormalization",
    .model = model,
    .bootstrap=bootstrap,
    .chunkSize=as.integer(500e3),
    .maxIter=as.integer(50),
    .acc=0.005
  )
})


setMethodS3("getAsteriskTags", "BasePositionNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Add model tag?
  model <- this$.model;
  if (model != "smooth.spline") {
    tags <- c(tags, model);
  }

  # Add bootstrap tag?
  if (this$.bootstrap) {
    bootstrapTag <- "B";
    tags <- c(tags, bootstrapTag);
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)



setMethodS3("getParameters", "BasePositionNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  params <- c(params, list(
    model = this$.model,
    bootstrap = this$.bootstrap,
    chunkSize = this$.chunkSize,
    maxIter = this$.maxIter,
    acc = this$.acc
  ));

  params;
}, private=TRUE)



setMethodS3("getDesignMatrix", "BasePositionNormalization", function(this, cells=NULL, model=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving design matrix");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  verbose && cat(verbose, "Model: ", model);

  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));

  key <- list(method="getDesignMatrix", class=class(this[1]), cells=cells, model=model, acs=list(fullname=getFullName(acs), checksum=getChecksum(acs)));
  dirs <- c("aroma.affymetrix", getChipType(acs));
  if (!force) {
    X <- loadCache(key=key, dirs=dirs);
    if (!is.null(X)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(X);
    }
  }

  verbose && enter(verbose, "Reading probe sequences");
  seqs <- readSequenceMatrix(acs, cells=cells, what="raw", 
                                                verbose=less(verbose, 5));
  verbose && cat(verbose, "Probe-sequence matrix:");
  verbose && str(verbose, seqs);
  verbose && exit(verbose);

  verbose && enter(verbose, "Building probe-position design matrix");
  X <- getProbePositionEffectDesignMatrix(seqs, verbose=less(verbose, 5));
  rm(seqs);

  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, X);
  verbose && cat(verbose, "RAM: ", object.size(X));
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Cache results?
  if (cache) {
    saveCache(X, key=key, dirs=dirs);
  }

  X;
}, private=TRUE)



setMethodS3("fitOne", "BasePositionNormalization", function(this, df, cells=NULL, ..., verbose=FALSE) {
  fitSubset <- function(df, cells=NULL, ..., verbose) {
    verbose && enter(verbose, "Retrieving probe sequences");

    verbose && enter(verbose, "Excluding cells with unknown probe sequences");
    # Exclude missing sequences
    n <- length(cells);
    acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 20));
    isMissing <- isMissing(acs, verbose=less(verbose, 10))[cells];
    cells <- cells[!isMissing];
    rm(isMissing);
    n2 <- length(cells);
    verbose && printf(verbose, "Removed %d (%.4f%%) missing sequences out of %d: ", n-n2, 100*(n-n2)/n, n);
    verbose && exit(verbose);

    verbose && enter(verbose, "Getting design matrix");
    X <- getDesignMatrix(this, cells=cells, verbose=less(verbose, 5));
    verbose && cat(verbose, "Design matrix:");
    verbose && str(verbose, X);
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading signals to fit");
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    y <- extractMatrix(df, cells=cells, drop=TRUE, verbose=less(verbose, 10));
    rm(cells);
    verbose && exit(verbose);

    if (shift != 0) {
      verbose && enter(verbose, "Shifting signals");
      verbose && cat(verbose, "Shift: ", shift);
      y <- y + shift;
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Log2 transforming signals");
    y <- log2(y);
    verbose && cat(verbose, "Target log2 probe signals:");
    verbose && str(verbose, y);
    verbose && exit(verbose);

    verbose && enter(verbose, "Keeping only finite data points");
    n <- length(y);
    # Fit only finite subset
    keep <- whichVector(is.finite(y));
    y <- y[keep];
    X$X <- X$X[keep,,drop=FALSE];
    rm(keep);
    n2 <- length(y);
    verbose && printf(verbose, "Removed %d (%.4f%%) non-finite data points out of %d: ", n-n2, 100*(n-n2)/n, n);
    verbose && exit(verbose);

    verbose && enter(verbose, "Fitting base-count model");
    verbose && cat(verbose, "Log2 signals:");
    verbose && str(verbose, y);
    verbose && cat(verbose, "Design matrix:");
    verbose && str(verbose, X);
    gc <- gc();
    verbose && print(verbose, gc);
    fit <- fitProbePositionEffects(y=y, seqs=X, verbose=less(verbose, 5));
    rm(y, X);
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  
    fit;
  } # fitSubset()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  if (!inherits(df, "AffymetrixCelFile")) {
    throw("Argument 'df' is not an AffymetrixCelFile: ", class(df)[1]);
  }

  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting normalization function for one array");
  verbose && cat(verbose, "Full name: ", getFullName(df));

  # Get algorithm parameters
  params <- getParameters(this, expand=FALSE, verbose=less(verbose, 5));
  model <- params$model;
  shift <- params$shift;
  bootstrap <- params$bootstrap;
  chunkSize <- params$chunkSize;
  maxIter <- params$maxIter;
  acc <- params$acc;

  # Is bootstrapping necessary?
  if (chunkSize >= length(cells)) {
    verbose && cat(verbose, "Bootstrapping not really needed.");
    bootstrap <- FALSE;
  }

  # Bootstrap?
  if (bootstrap) {
    verbose && enter(verbose, "Fitting model using bootstrap");
    verbose && cat(verbose, "Max number of iterations: ", maxIter);
    verbose && cat(verbose, "Accuracy threshold: ", acc);
throw("Yet not implemented!");

    bb <- 0;
    fit <- NULL;
    deltaMax <- Inf;
    while (deltaMax > acc && bb < maxIter) {
      bb <- bb + 1;
      verbose && enter(verbose, sprintf("Bootstrap iteration #%d", as.integer(bb)));

      verbose && enter(verbose, "Fitting model to subset of data");
      verbose && cat(verbose, "Chunk size: ", chunkSize);
      subset <- sample(1:length(cells), size=chunkSize);
      cellsChunk <- cells[subset];
      rm(subset);
      cellsChunk <- sort(cellsChunk);
      verbose && cat(verbose, "Cells:");
      verbose && str(verbose, cellsChunk);
      fitB <- fitSubset(df, cells=cellsChunk, verbose=verbose);
      rm(cellsChunk);
      verbose && exit(verbose);

      verbose && enter(verbose, "Updating estimates");
      if (is.null(fit)) {
        fit <- fitB;
      } else {
        # Update parameters for each subfit
        deltaMax <- 0;
      }
      rm(fitB);
      verbose && exit(verbose);

      verbose && printf(verbose, "deltaMax < threshold: %.6f < %.6f\n", deltaMax, acc);
      verbose && printf(verbose, "iteration < maxIter: %d < %d\n", as.integer(bb), as.integer(maxIter));

      verbose && exit(verbose);
    } # while()
    converged <- (deltaMax <= acc);

    verbose && enter(verbose, "Creating final fit");
    verbose && exit(verbose);

    fit$bootstrap <- list(iter=as.integer(bb), maxIter=maxIter, converged=converged);

    verbose && exit(verbose);
  } else {
    fit <- fitSubset(df, cells=cells, verbose=verbose);
  } # if (bootstrap)

  verbose && exit(verbose);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "BasePositionNormalization", function(this, fit, cells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Predicting model for one array");

  # Get algorithm parameters
  params <- getParameters(this, expand=FALSE, verbose=less(verbose, 5));
  model <- params$model;

  verbose && enter(verbose, "Retrieving probe sequences");
  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  seqs <- readSequenceMatrix(acs, cells=cells, what="raw", verbose=less(verbose, 5));
  rm(acs);
  verbose && cat(verbose, "Probe-sequence matrix:");
  verbose && str(verbose, seqs);
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  verbose && enter(verbose, "Predicting mean log2 probe signals");
  mu <- predict(fit, seqs=seqs, verbose=less(verbose, 5));
  rm(seqs);
  verbose && str(verbose, "mu:");
  verbose && str(verbose, mu);
  verbose && exit(verbose);

  # Sanity check
  if (length(mu) != length(cells)) {
    throw("Internal error. Number of estimated means does not match the number of cells to be updated: ", length(mu), " != ", length(cells));
  }

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  mu;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-07-21
# o BENCHMARKING: For a GenomeWideSNP_6,Full, the BPN peaks at 5.9GB RAM.
#   This happens while fitting the model.  Prediction peaks at 3.2GB RAM.
# o Now getDesignMatrix() caches results to file.
# o Created from BaseCountNormalization.R.
############################################################################
