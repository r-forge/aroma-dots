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
#   \item{df}{The degrees of freedom of the model.}
#   \item{.fitMethod}{...}
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
setConstructorS3("BasePositionNormalization", function(..., model=c("smooth.spline"), df=5, .fitMethod=c("lm.fit", "solve")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'df':
  df <- Arguments$getInteger(df, range=c(1,1e3));

  # Argument '.fitMethod':
  .fitMethod <- match.arg(.fitMethod);


  extend(AbstractProbeSequenceNormalization(...), "BasePositionNormalization",
    .model = model,
    .df = df,
    .fitMethod = .fitMethod
  )
})


setMethodS3("getAsteriskTags", "BasePositionNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Add model tag?
  model <- this$.model;
  if (model != "smooth.spline") {
    tags <- c(tags, model);
  }

  # Add df tag?
  df <- this$.df;
  if (df != 5) {
    tags <- c(tags, sprintf("df=%d", df));
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
    df = this$.df,
    fitMethod = this$.fitMethod
  ));

  params;
}, private=TRUE)



setMethodS3("getDesignMatrix", "BasePositionNormalization", function(this, cells=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
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

  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=FALSE, verbose=less(verbose, 1));
  model <- params$model;
  df <- params$df;
  verbose && cat(verbose, "Model: ", model);
  verbose && cat(verbose, "Degrees of freedom: ", df);
  verbose && exit(verbose);

  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));

  key <- list(
    method="getDesignMatrix", class=class(this[1]), 
    cells=cells, 
    model=model, df=df, 
    acs=list(fullname=getFullName(acs), checksum=getChecksum(acs))
  );
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
  verbose && cat(verbose, "Degrees of freedom: ", df);
  X <- getProbePositionEffectDesignMatrix(seqs, df=df, verbose=less(verbose, 5));
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



setMethodS3("getNormalEquations", "BasePositionNormalization", function(this, df=NULL, cells=NULL, transform=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  if (!inherits(df, "AffymetrixCelFile")) {
    throw("Argument 'df' is not an AffymetrixCelFile: ", class(df)[1]);
  }

  # Argument 'transform':
  if (!is.null(transform)) {
    if (!is.function(transform)) {
      throw("Argument 'transform' is not a function: ", mode(transform)[1]);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting annotation data files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving cell sequence annotation data file");
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 20));
  verbose && exit(verbose);

  # Expand 'cells'?
  if (is.null(cells)) {
    cells <- 1:nbrOfCells(acs);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying subset of cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying subset of cells that can be fitted");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  n0 <- length(cells);

  verbose && enter(verbose, "Excluding cells with unknown probe sequences");
  isMissing <- isMissing(acs, verbose=less(verbose, 10))[cells];
  cells <- cells[!isMissing];
  rm(isMissing);
  n1 <- length(cells);
  verbose && printf(verbose, "Removed %d (%.2f%%) missing sequences out of %d\n", n0-n1, 100*(n0-n1)/n0, n0);
  verbose && exit(verbose);


  verbose && enter(verbose, "Identifying cells with missing data");
  verbose && enter(verbose, "Reading signals");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  y <- extractMatrix(df, cells=cells, drop=TRUE, verbose=less(verbose, 10));

  if (!is.null(transform)) {
    verbose && enter(verbose, "Transforming signals");
    verbose && cat(verbose, "Signals before transformation:");
    verbose && str(verbose, y);
    y <- transform(y);
    verbose && exit(verbose);
  }
  verbose && cat(verbose, "Signals to be fitted:");
  verbose && str(verbose, y);
  verbose && exit(verbose);
    
  verbose && enter(verbose, "Excluding non-finite data points");
  # Fit only finite subset
  keep <- whichVector(is.finite(y));
  y <- y[keep];
  cells <- cells[keep];
  rm(keep);
  n2 <- length(cells);
  verbose && printf(verbose, "Removed %d (%.2f%%) non-finite data points out of %d\n", n1-n2, 100*(n1-n2)/n1, n1);
  verbose && exit(verbose);

  verbose && printf(verbose, "Removed in total %d (%.2f%%) cells out of %d\n", n0-n2, 100*(n0-n2)/n0, n0);
  verbose && str(verbose, cells);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting design matrix for subset
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting design matrix");
  res <- getDesignMatrix(this, cells=cells, verbose=less(verbose, 5));
  X <- res$X;

##  XX1 <<- X; yy1 <<- y;

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

  gc <- gc();
  verbose && print(verbose, gc);

  # Sanity check
  stopifnot(nrow(X) == length(cells));
  stopifnot(nrow(X) == length(y));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculating cross products X'X and X'y
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating cross product X'X");
  xtx <- crossprod(X);  
  verbose && str(verbose, xtx);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating cross product X'y");
  xty <- crossprod(X, y);
  verbose && str(verbose, xty);
  verbose && exit(verbose);

  rm(X);

  res <- list(xtx=xtx, xty=xty, n0=n0, n1=n1, n2=n2, cells=cells, map=map, B=B, factors=factors, X=X, y=y);
  rm(xtx, xty, cells);

  res;
}, protected=TRUE)




setMethodS3("fitOne", "BasePositionNormalization", function(this, df, ..., verbose=FALSE) {
  fitSubset <- function(df, cells=NULL, shift=0, ..., verbose) {
    verbose && enter(verbose, "Retrieving probe sequences");

    verbose && enter(verbose, "Excluding cells with unknown probe sequences");
    # Exclude missing sequences
    n <- length(cells);
    acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 20));
    isMissing <- isMissing(acs, verbose=less(verbose, 10))[cells];
    cells <- cells[!isMissing];
    rm(isMissing);
    n2 <- length(cells);
    verbose && printf(verbose, "Removed %d (%.2f%%) missing sequences out of %d\n", n-n2, 100*(n-n2)/n, n);
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
    verbose && cat(verbose, "Log2 probe signals:");
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
    verbose && printf(verbose, "Removed %d (%.2f%%) non-finite data points out of %d\n", n-n2, 100*(n-n2)/n, n);
    verbose && exit(verbose);

    verbose && enter(verbose, "Fitting base-count model");
    verbose && cat(verbose, "Log2 signals:");
    verbose && str(verbose, y);
    verbose && cat(verbose, "Design matrix:");
    verbose && str(verbose, X);
##  XX2 <<- X$X; yy2 <<- y;

    gc <- gc();
    verbose && print(verbose, gc);
    fit <- fitProbePositionEffects(y=y, seqs=X, verbose=less(verbose, 5));
    rm(y, X);
    fit$algorithm <- "lm.fit";
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  
    fit;
  } # fitSubset()

  fitSubset2 <- function(df, cells=NULL, shift=0, ..., verbose) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Local functions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    transform <- function(y, ...) {
      if (shift != 0) {
        verbose && enter(verbose, "Shifting signals");
        verbose && cat(verbose, "Shift: ", shift);
        y <- y + shift;
        verbose && exit(verbose);
      }
      verbose && enter(verbose, "Log2 transforming signals");
      y <- log2(y);
      verbose && exit(verbose);
      y;
    } # transform()

    ne <- getNormalEquations(this, df=df, cells=cells, transform=transform, verbose=verbose);
    verbose && cat(verbose, "Normal equations:");
    verbose && str(verbose, ne);

print(sum(ne$X)); print(sum(ne$y));

    map <- ne$map;
    B <- ne$B;
    factors <- ne$factors;

    xtx <- ne$xtx;
    xty <- ne$xty;
    rm(ne);

    coefs <- solve(xtx, xty);
print(coefs);
    coefs <- as.vector(coefs);
    verbose && cat(verbose, "Coeffients:")
    verbose && print(verbose, coefs);

    params <- list();
    intercept <- TRUE;
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
    fit <- list(params=params, map=map, B=B, algorithm="solve");
    class(fit) <- "ProbePositionEffects";

#    params <- list(coefs=coefs);
#    fit <- list(params=params, map=map, B=B);
#    class(fit) <- c("fitSubset2", class(fit));
    verbose && str(verbose, fit);

    fit;
  } # fitSubset2()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  if (!inherits(df, "AffymetrixCelFile")) {
    throw("Argument 'df' is not an AffymetrixCelFile: ", class(df)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting normalization function for one array");
  verbose && cat(verbose, "Full name: ", getFullName(df));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));
  gc <- gc();
  verbose && print(verbose, gc);

  units <- params$unitsToFit;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  cells <- params$cellsToFit;
  stratifyBy <- params$typesToFit;
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  shift <- params$shift;
  verbose && cat(verbose, "Shift: ", shift);

  # Model parameters (only for display; retrieve elsewhere)
  verbose && cat(verbose, "Model: ", params$model);
  verbose && cat(verbose, "Degrees of freedom: ", params$df);

  # Method for fitting linear model
  fitMethod <- params$fitMethod;
  verbose && cat(verbose, "Algorithm for fitting model: ", fitMethod);
  verbose && exit(verbose);

  if (fitMethod == "lm.fit") {
    verbose && enter(verbose, "Exact fitting of model using lm.fit");
    verbose && cat(verbose, "NOTE: This approach is not bounded in memory. For larger chip types it may peak at 4-6GB of RAM.");
    fit <- fitSubset(df, cells=cells, shift=shift, verbose=verbose);
    verbose && exit(verbose);
  } else if (fitMethod == "solve") {
    verbose && enter(verbose, "Exact fitting of model by incrementally building the normal equations (X'X = X'y) and then solve it");
    fit <- fitSubset2(df, cells=cells, shift=shift, verbose=verbose);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "BasePositionNormalization", function(this, fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Predicting model for one array");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));

  units <- params$unitsToUpdate;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  cells <- params$cellsToUpdate;
  stratifyBy <- params$typesToUpdate;
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  shift <- params$shift;
  verbose && cat(verbose, "Shift: ", shift);

  # Other model parameters
  model <- params$model;
  verbose && cat(verbose, "Model: ", model);
  verbose && exit(verbose);

  nbrOfCells <- length(cells);


  verbose && enter(verbose, "Retrieving probe sequences");
  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  seqs <- readSequenceMatrix(acs, cells=cells, what="raw", verbose=less(verbose, 5));
  rm(acs, cells);
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
  verbose && summary(verbose, mu);
  verbose && exit(verbose);

  # Sanity check
  if (length(mu) != nbrOfCells) {
    throw("Internal error. Number of estimated means does not match the number of cells to be updated: ", length(mu), " != ", nbrOfCells);
  }

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  mu;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-11-29
# o Added first step toward supporting fitting the linear model in
#   bounded memory.  This is done by setting up the normal equations and
#   using solve(xtx, xty) to estimate the parameters.
#   TODO: Build up the NE incrementally.
# o Dropped the bootstrapping framework.
# 2008-07-29
# o Added support for specifying the degrees of freedom ('df') of the model.
# 2008-07-28
# o Updated to work with newer ProbeLevelTransform3.
# 2008-07-21
# o BENCHMARKING: For a GenomeWideSNP_6,Full, the BPN peaks at 5.9GB RAM.
#   This happens while fitting the model.  Prediction peaks at 3.2GB RAM.
# o Now getDesignMatrix() caches results to file.
# o Created from BaseCountNormalization.R.
############################################################################
