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
setConstructorS3("BaseCountNormalization", function(dataSet=NULL, ..., typesToUpdate="pm", subsetToUpdate=NULL, model=c("robustSmoothSpline", "lm"), typesToFit=typesToUpdate, subsetToFit="-XY") {
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

    cdf <- getCdf(dataSet);

    # Argument 'typesToUpdate':
    if (!is.null(typesToUpdate)) {
      typesToUpdate <- match.arg(typesToUpdate, choices=c("pm", "mm", "pmmm"));
    }

    # Argument 'subsetToUpdate':
    if (is.null(subsetToUpdate)) {
    } else if (is.character(subsetToUpdate)) {
      if (subsetToUpdate %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToUpdate': ", subsetToUpdate);
      }
      extraTags <- c(extraTags, subsetToUpdate=subsetToUpdate);
    } else {
      nbrOfCells <- nbrOfCells(cdf);
      subsetToUpdate <- Arguments$getIndices(subsetToUpdate, range=c(1, nbrOfCells));
      subsetToUpdate <- unique(subsetToUpdate);
      subsetToUpdate <- sort(subsetToUpdate);
    }

    # Argument 'model':
    model <- match.arg(model);

    # Argument 'typesToFit':
    if (!is.null(typesToFit)) {
      typesToFit <- match.arg(typesToFit, choices=c("pm", "mm", "pmmm"));
    }

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
    .typesToUpdate = typesToUpdate,
    .subsetToUpdate = subsetToUpdate,
    .model = model,
    .typesToFit = typesToFit,
    .subsetToFit = subsetToFit,
    .extraTags = extraTags
  )
})


setMethodS3("clearCache", "BaseCountNormalization", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".subsetToUpdateExpanded", ".subsetToFitExpanded")) {
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




setMethodS3("getSubsetTo", "BaseCountNormalization", function(this, what=c("fit", "update"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);

  field <- sprintf(".subsetTo%s", capitalize(what));
  fieldExpanded <- paste(field, "Expanded", sep="");
  typesField <- sprintf(".typesTo%s", capitalize(what));

  subset <- this[[field]];
  stratifyBy <- this[[typesField]];

  # Expand?
  if (is.character(subset)) {
    cells <- this[[fieldExpanded]];
    if (is.null(cells)) {
      dataSet <- getInputDataSet(this);
      cdf <- getCdf(dataSet);
      units <- subset;
      cells <- getSubsetOfCellIndices(cdf, units=units, stratifyBy=stratifyBy, ...);
      this[[fieldExpanded]] <- cells;
    }
  } else if (is.numeric(subset)) {
    cells <- subset;
  } else if (is.null(subset)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    if (is.null(stratifyBy)) {
      cells <- seq(length=nbrOfCells(cdf));
    } else {
      cells <- getSubsetOfCellIndices(cdf, stratifyBy=stratifyBy, ...);
    }
    this[[fieldExpanded]] <- cells;
  } else {
    throw("Internal error. This statment should never be reached: ", mode(subset));
  }

  cells;
}, private=TRUE);


setMethodS3("getSubsetToUpdate", "BaseCountNormalization", function(this, ...) {
  getSubsetTo(this, what="update", ...);
}, protected=TRUE);


setMethodS3("getSubsetToFit", "BaseCountNormalization", function(this, ...) {
  getSubsetTo(this, what="fit", ...);
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
    typesToUpdate = this$.typesToUpdate,
    subsetToUpdate = this$.subsetToUpdate,
    model = this$.model,
    typesToFit = this$.typesToFit,
    subsetToFit = this$.subsetToFit
  ));

  # Expand?
  if (expand) {
    params$subsetToUpdate <- getSubsetToUpdate(this);
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

  # Get subset to update
  subsetToUpdate <- params$subsetToUpdate;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- nbrOfArrays(ds);
  nbrOfCells <- nbrOfCells(cdf);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);

  designMatrix <- NULL;
  paramsShort <- NULL;
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
      # Generate design matrix for all probes of interest?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(designMatrix)) {
        verbose && enter(verbose, "Count nucleotide bases for *all* cells");
        verbose && cat(verbose, "Chip type: ", getChipType(cdf));
        designMatrix <- countBases(this, verbose=less(verbose, 5));
        verbose && cat(verbose, "Nucleotide base counts:");
        verbose && str(verbose, designMatrix);

        designMatrix[,1] <- as.integer(1);
        verbose && cat(verbose, "Design matrix:");
        verbose && str(verbose, designMatrix);

        # Identify missing sequences (to be excluded from fit and updates)
        missingSeqs <- which(is.na(designMatrix[,4]));

        # Update subset of cell indices for fitting and updating
        subsetToFit <- setdiff(subsetToFit, missingSeqs);
        subsetToUpdate <- setdiff(subsetToUpdate, missingSeqs);

        verbose && cat(verbose, "Cell indices used for fitting:");
        verbose && str(verbose, subsetToFit);
        verbose && cat(verbose, "Cell indices to be updated:");
        verbose && str(verbose, subsetToUpdate);

        rm(missingSeqs); # Not needed anymore
        gc <- gc();
        verbose && exit(verbose);
      }

      if (is.null(paramsShort)) {
        # Precalculate some model fit parameters
        verbose && enter(verbose, "Compressing model parameter to a short format");
        paramsShort <- params;
        paramsShort$subsetToFit <- NULL;
        paramsShort$subsetToUpdate <- NULL;
        paramsShort$subsetToFitIntervals <- seqToIntervals(subsetToFit);
        paramsShort$subsetToUpdateIntervals <- seqToIntervals(subsetToUpdate);
        verbose && exit(verbose);
      }


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Setting up model fit parameters
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      modelFit <- list(
        paramsShort=paramsShort
      );



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase 0: Fit base-count effect for target?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(muT)) {
        verbose && enter(verbose, "Estimating base-count effects for target");

        verbose && enter(verbose, "Getting target signals to fit the model");
        dfT <- getTargetFile(this, verbose=less(verbose, 5));
        yT <- extractMatrix(dfT, cells=subsetToFit, drop=TRUE);
        rm(dfT);

        verbose && cat(verbose, "Target log2 probe signals:");
        yT <- log2(yT);
        verbose && str(verbose, yT);
        gc <- gc();
        verbose && exit(verbose);

        # Fit only finite subset
        keep <- which(is.finite(yT));
        yT <- yT[keep];
        subsetToFitT <- subsetToFit[keep];
        rm(keep);
        gc <- gc();

        verbose && enter(verbose, "Fitting base-count model");
        X <- designMatrix[subsetToFitT,,drop=FALSE];
        rm(subsetToFitT);
        verbose && cat(verbose, "Design matrix:");
        verbose && str(verbose, X);
        gc <- gc();
        verbose && print(verbose, gc);
        fitT <- fitBaseCounts(yT, X=X, model=params$model, verbose=less(verbose, 5));
        rm(yT, X);
        verbose && print(verbose, fitT);
        verbose && exit(verbose);


        verbose && enter(verbose, "Target mean log2 probe signals:");
        X <- designMatrix[subsetToUpdate,,drop=FALSE];
        verbose && cat(verbose, "Design matrix:");
        verbose && str(verbose, X);
        muT <- predictBaseCounts(fitT, X=X, model=params$model);
        rm(fitT, X);
        verbose && str(verbose, "muT:");
        verbose && str(verbose, muT);
        if (length(muT) != length(subsetToUpdate)) {
          throw("Internal error. Number of estimated means does not match the number of cells to be updated: ", length(mu), " != ", length(subsetToUpdate));
        }
        verbose && exit(verbose);

        gc <- gc();

        verbose && exit(verbose);
      } # if (is.null(muT))



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase I: Fit base-count effect for the current array
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Getting signals used to fit the model");
      y <- extractMatrix(df, cells=subsetToFit, drop=TRUE);
      y <- log2(y);
      verbose && cat(verbose, "Log2 probe signals:");
      verbose && str(verbose, y);
      gc <- gc();
      verbose && exit(verbose);

      # Find finite subset
      keep <- which(is.finite(y));
      y <- y[keep];
      subsetToFitKK <- subsetToFit[keep];
      rm(keep);
      gc <- gc();

      verbose && enter(verbose, "Fitting base-count model");
      X <- designMatrix[subsetToFitKK,,drop=FALSE];
      verbose && cat(verbose, "Design matrix:");
      verbose && str(verbose, X);
      rm(subsetToFitKK);
      gc <- gc();
      verbose && print(verbose, gc);
      fit <- fitBaseCounts(y, X=X, model=params$model,
                                                 verbose=less(verbose, 5));
      rm(y, X);
      verbose && print(verbose, fit);
      modelFit$fit <- fit;
      verbose && exit(verbose);

      # Store model fit 
      verbose && enter(verbose, "Saving model fit");
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename, 
                                                      path=outputPath, ...);
      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      rm(modelFit);
      verbose && exit(verbose);

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase II: Normalize current array toward target
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading probe signals:");
      y <- extractMatrix(df, cells=subsetToUpdate, drop=TRUE);
      verbose && str(verbose, y);
      verbose && summary(verbose, y);
      verbose && exit(verbose);

      # Find finite subset for this array
      keep <- which(is.finite(y));
      subsetToUpdateKK <- subsetToUpdate[keep];
      y <- y[keep];
      gc <- gc();

      X <- designMatrix[subsetToUpdateKK,,drop=FALSE];
      verbose && cat(verbose, "Design matrix:");
      verbose && str(verbose, X);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && enter(verbose, "Predicting mean log2 probe signals:");
      mu <- predictBaseCounts(fit, X=X, model=params$model);
      rm(fit, X);
      verbose && str(verbose, "mu:");
      verbose && str(verbose, mu);
      if (length(mu) != length(y)) {
        throw("Internal error. Number of estimated means does not match the number of data points: ", length(mu), " != ", length(y));
      }
      verbose && exit(verbose);


      verbose && enter(verbose, "Discrepancy scale factors:");
      rho <- (muT[keep]-mu);
      rm(mu, keep);
      summary(verbose, rho);
      rho <- 2^rho;
      summary(verbose, rho);

      # Update only subset with "finite" corrections
      keep <- which(is.finite(rho));
      rho <- rho[keep];
      subsetToUpdateKK <- subsetToUpdateKK[keep];
      y <- y[keep];
      rm(keep);
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Normalizing probe signals:");
      y <- rho * y;
      rm(rho);
      verbose && str(verbose, y);
      verbose && summary(verbose, y);
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
      updateCel(pathname, indices=subsetToUpdateKK, intensities=y, verbose=verbose2);
      rm(y, subsetToUpdateKK, verbose2);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
    }

    # Validating by retrieving calibrated data file
    dfC <- newInstance(df, pathname);

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
# 2008-07-16
# o Added support for fitting and updating subsets of cells and types of
#   probes according to the CDF.
# 2008-07-11
# o Updated countBases() to use new AromaCellSequenceFile class.
# 2008-06-22
# o Created from AllelicCrosstalkCalibration.R.
############################################################################
