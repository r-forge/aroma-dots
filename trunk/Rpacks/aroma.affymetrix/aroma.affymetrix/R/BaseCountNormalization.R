###########################################################################/**
# @RdocClass BaseCountNormalization
#
# @title "The BaseCountNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in the number of
#  A, C, G, and T:s in the probe sequences.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{model}{A @character string specifying the model used to fit 
#     the base-count effects.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
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
setConstructorS3("BaseCountNormalization", function(dataSet=NULL, ..., model=c("robustSmoothSpline", "lm"), subsetToFit="-XY") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extraTags <- NULL;

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", 
                                                          class(dataSet)[1]);
    }

    # Argument 'model':
    model <- match.arg(model);

    cdf <- getCdf(dataSet);

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      if (subsetToFit %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToFit': ", subsetToFit);
      }
      extraTags <- c(extraTags, subsetToFit=subsetToFit);
    } else {
      nbrOfCells <- nbrOfCells(cdf);
      subsetToFit <- Arguments$getIndices(subsetToFit, range=c(1, nbrOfCells));
      subsetToFit <- unique(subsetToFit);
      subsetToFit <- sort(subsetToFit);
    }
  }


  extend(ProbeLevelTransform(dataSet=dataSet, ...), "BaseCountNormalization",
    .model = model,
    .subsetToFit = subsetToFit,
    .extraTags = extraTags
  )
})


setMethodS3("clearCache", "BaseCountNormalization", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".subsetToFitExpanded")) {
    this[[ff]] <- NULL;
  }

  # Then for this object 
  NextMethod("clearCache", object=this, ...);
})


setMethodS3("getAsteriskTags", "BaseCountNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)



setMethodS3("getSubsetToFit", "BaseCountNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  subsetToFit <- this$.subsetToFit;

  # Expand?
  if (is.character(subsetToFit)) {
    if (subsetToFit %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToFit: ", subsetToFit);

      # Look up in cache
      subset <- this$.subsetToFitExpanded;
      if (is.null(subset)) {
        dataSet <- getInputDataSet(this);
        cdf <- getCdf(dataSet);
  
        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);
  
        # Identify units to be excluded
        if (subsetToFit == "-X") {
          subset <- getUnitsOnChromosome(gi, 23, .checkArgs=FALSE);
        } else if (subsetToFit == "-Y") {
          subset <- getUnitsOnChromosome(gi, 24, .checkArgs=FALSE);
        } else if (subsetToFit == "-XY") {
          subset <- getUnitsOnChromosome(gi, 23:24, .checkArgs=FALSE);
        }
  
        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, subset);

        # Identify the cell indices for these units
        subset <- getCellIndices(cdf, units=subset, 
                                 useNames=FALSE, unlist=TRUE);
        verbose && cat(verbose, "Cells to exclude: ");
        verbose && str(verbose, subset);
  
        # The cells to keep
        subset <- setdiff(1:nbrOfCells(cdf), subset);
  
        verbose && cat(verbose, "Cells to include: ");
        verbose && str(verbose, subset);

        # Store
        this$.subsetToFitExpanded <- subset;
      }

      subsetToFit <- subset;
      rm(subset);

      verbose && exit(verbose);
    }
  }

  subsetToFit;
}, protected=TRUE);


setMethodS3("getTargetFile", "BaseCountNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  dataSet <- getInputDataSet(this);
  dfR <- getAverageFile(dataSet, verbose=less(verbose, 25));

  dfR;
})


setMethodS3("getAromaCellSequenceFile", "BaseCountNormalization", function(this, ..., force=FALSE) {
  aps <- this$.aps;

  if (is.null(aps)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf);
    aps <- AromaCellSequenceFile$byChipType(chipType, ...);
    this$.aps <- aps;
  }

  aps;
}, protected=TRUE)


setMethodS3("countBases", "BaseCountNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  counts <- this$.baseCounts;

  if (force || is.null(counts)) {
    aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 20));
    counts <- countBases(aps, verbose=less(verbose, 5));
    this$.baseCounts <- counts;
  }

  counts;
}, protected=TRUE);


setMethodS3("getParameters", "BaseCountNormalization", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, expand=expand, ...);

  params <- c(params, list(
    model = this$.model,
    subsetToFit = this$.subsetToFit
  ));

  # Expand?
  if (expand) {
    params$subsetToFit <- getSubsetToFit(this);
  }

  params;
}, private=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Calibrates the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already calibrated is re-calibrated, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "BaseCountNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'y' and 'X' must not contain NAs.
  fitBaseCounts <- function(y, X, subset=NULL, model=c("robustSmoothSpline", "lm"), ...) {
    # Argument 'y':

    # Argument 'X':
    if (nrow(X) != length(y)) {
      throw("The number of rows in design matrix 'X' does not match the number of observations in 'y': ", nrow(X), " != ", length(y));
    }

    if (!is.null(subset)) {
      y <- y[subset];
      X <- X[subset,,drop=FALSE];
      gc <- gc();
    }

    # Argument 'model':
    model <- match.arg(model);
print(model);

    if (model == "lm") {
      require("stats") || throw("Package not loaded: stats");
      fitFcn <- function(X, y, ...) {
        fit <- stats::lm.fit(x=X, y=y, ...);
        # Remove redundant parameters
        for (ff in c("residuals", "effects", "fitted.values", "qr")) {
          fit[[ff]] <- NULL;
        }
        fit;
      }
    } else if (model == "robustSmoothSpline") {
      require("aroma.light") || throw("Package not loaded: aroma.light");
      fitFcn <- function(X, y, ...) {
        fits <- list();
        for (cc in 1:ncol(X)) {
          # Fit effect of term #cc
          if (cc == 1) {
            mu <- median(y);
            fit <- list(mu=mu);
          } else {
            fit <- aroma.light::robustSmoothSpline(x=X[,cc], y=y, ...);
            # Remove redundant parameters
            for (ff in c("x", "y", "w", "yin", "lev")) {
              fit[[ff]] <- NULL;
            }
            mu <- predict(fit, x=X[,cc])$y;
          }

          # Remove the effect of term #cc
          y <- y - mu;
          rm(mu);

          fits[[cc]] <- fit;
          rm(fit);
        }
        fits;
      }
    }

    fit <- fitFcn(X, y);

    fit;
  } # fitBaseCounts()


  predictBaseCounts <- function(fit, X, model=c("robustSmoothSpline", "lm"), ...) {
    # Argument 'model':
    model <- match.arg(model);
print(model);

    if (model == "lm") {
      predictFcn <- function(fit, X, ...) {
        coefs <- coefficients(fit);
        coefs <- as.matrix(coefs);
        yPred <- X %*% coefs;
        yPred;
      } # predictFcn()
    } else if (model == "robustSmoothSpline") {
      predictFcn <- function(fit, X, ...) {
        fits <- fit;
        yPred <- double(nrow(X));
        for (cc in 1:ncol(X)) {
          fit <- fits[[cc]];
          if (cc == 1) {
            mu <- fit$mu;
          } else {
            idxs <- which(is.finite(X[,cc]));
            mu <- double(nrow(X));
            mu[idxs] <- predict(fit, x=X[idxs,cc])$y;
            rm(idxs);
          }
          str(mu);
          yPred <- yPred + mu;
          rm(fit, mu);
        } # for (cc ...)
        yPred;
      } # predictFcn()
    }

    yPred <- predictFcn(fit, X);
    
    yPred;
  } # predictBaseCounts()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalization data set for probe sequence base count effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Get subset to fit
  subsetToFit <- params$subsetToFit;

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "process.BaseCountNormalization";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate some model fit parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Compressing model parameter to a short format");
  paramsShort <- params;
  paramsShort$subsetToFit <- NULL;
#  paramsShort$subsetToFitIntervals <- seqToIntervals(params$subsetToFit);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- nbrOfArrays(ds);
  nbrOfCells <- nbrOfCells(cdf);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);

  designMatrix <- NULL;
  muT <- NULL;
  hasSeq <- NULL;
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already calibrated?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ", pathname);
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Calibrating
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      callHooks(sprintf("%s.onBegin", hookName), df=df, ...);

      modelFit <- list(
        paramsShort=paramsShort
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Count nucleotide bases for this chip type?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(designMatrix)) {
        verbose && enter(verbose, "Count nucleotide bases");
        verbose && cat(verbose, "Chip type: ", getChipType(cdf));
        designMatrix <- countBases(this, verbose=less(verbose, 5));
        verbose && cat(verbose, "Nucleotide base counts:");
        verbose && str(verbose, designMatrix);

        # Identify subset of cells to be used for fitting
        hasSeq <- which(!is.na(designMatrix[,4]));
        subsetToFit <- intersect(subsetToFit, hasSeq);
        rm(hasSeq);
        verbose && cat(verbose, "Subset of cells to be fitted:");
        verbose && str(verbose, subsetToFit);

        designMatrix[,1] <- as.integer(1);
        verbose && cat(verbose, "Design matrix:");
        verbose && str(verbose, designMatrix);
        gc <- gc();
        verbose && exit(verbose);
      }

      if (is.null(muT)) {
        verbose && enter(verbose, "Estimating base-count effects for target");

        verbose && enter(verbose, "Getting target signals to fit the model");
        dfT <- getTargetFile(this, verbose=less(verbose, 5));
        yT <- readRawData(dfT, indices=subsetToFit, fields="intensities", 
                                                                  drop=TRUE);
        rm(dfT);

        verbose && cat(verbose, "Target log2 probe signals:");
        yT <- log2(yT);
        verbose && str(verbose, yT);
        gc <- gc();
        verbose && exit(verbose);

        # Fit only finite subset
        subset <- which(is.finite(yT));
        yT <- yT[subset];
        gc <- gc();

        verbose && enter(verbose, "Fitting base-count model");
        subset <- subsetToFit[subset];
        X <- designMatrix[subset,,drop=FALSE];
        rm(subset);
        verbose && cat(verbose, "Design matrix:");
        verbose && str(verbose, X);
        gc <- gc();
        verbose && print(verbose, gc);
        fitT <- fitBaseCounts(yT, X=X, model=params$model, verbose=less(verbose, 5));
        rm(yT, X);
        verbose && print(verbose, fitT);
        verbose && exit(verbose);

        verbose && enter(verbose, "Target mean log2 probe signals:");
        muT <- predictBaseCounts(fitT, X=designMatrix, model=params$model);
        rm(fitT);
        verbose && str(verbose, "muT:");
        verbose && str(verbose, muT);
        if (length(muT) != nbrOfCells) {
          throw("Internal error. Number of estimated means does not match the number of cells on the array: ", length(mu), " != ", nbrOfCells);
        }
        verbose && exit(verbose);

        gc <- gc();

        verbose && exit(verbose);
      } # if (is.null(muT))


      verbose && enter(verbose, "Getting signals used to fit the model");
      y <- readRawData(df, indices=subsetToFit, fields="intensities", 
                                                                  drop=TRUE);
      y <- log2(y);
      verbose && cat(verbose, "Log2 probe signals:");
      verbose && str(verbose, y);
      gc <- gc();
      verbose && exit(verbose);

      # Find finite subset
      subset <- which(is.finite(y));
      y <- y[subset];
      gc <- gc();

      verbose && enter(verbose, "Fitting base-count model");
      subset <- subsetToFit[subset];
      X <- designMatrix[subset,,drop=FALSE];
      verbose && cat(verbose, "Design matrix:");
      verbose && str(verbose, X);
      rm(subset);
      gc <- gc();
      verbose && print(verbose, gc);
      fit <- fitBaseCounts(y, X=X, model=params$model, verbose=less(verbose, 5));
      rm(y, X);
      verbose && print(verbose, fit);
      modelFit$fit <- fit;
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store model fit 
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Saving model fit");
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      rm(modelFit);
      verbose && exit(verbose);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Normalize data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Mean log2 probe signals:");
      mu <- predictBaseCounts(fit, X=designMatrix, model=params$model);
      rm(fit);
      verbose && str(verbose, "mu:");
      verbose && str(verbose, mu);
      if (length(mu) != nbrOfCells) {
        throw("Internal error. Number of estimated means does not match the number of cells on the array: ", length(mu), " != ", nbrOfCells);
      }
      verbose && exit(verbose);


      verbose && enter(verbose, "Discrepancy scale factors:");
      rho <- (muT-mu);
      rm(mu);
      summary(verbose, rho);
      rho <- 2^rho;
      summary(verbose, rho);

      # Update only "finite" subset
      subset <- which(is.finite(rho));
      rho <- rho[subset];
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Reading probe signals:");
      yAll <- readRawData(df, fields="intensities", drop=TRUE);
      verbose && str(verbose, yAll);
      verbose && summary(verbose, yAll);
      verbose && exit(verbose);

      verbose && enter(verbose, "Normalizing probe signals:");
      yAll[subset] <- rho * yAll[subset];
      rm(rho, subset);
      verbose && str(verbose, yAll);
      verbose && summary(verbose, yAll);
      verbose && exit(verbose);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing normalized data");
    
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      updateCel(pathname, intensities=yAll, verbose=verbose2);
      rm(yAll, verbose2);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
    }

    # Validating by retrieving calibrated data file
    dfC <- newInstance(df, pathname);

#    callHooks(sprintf("%s.onExit", hookName), df=df, dfC=dfC, ...);

    rm(df, dfC);

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  # Garbage collect
  rm(ds);
  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})



############################################################################
# HISTORY:
# 2008-07-11
# o Updated countBases() to use new AromaCellSequenceFile class.
# 2008-06-22
# o Created from AllelicCrosstalkCalibration.R.
############################################################################
