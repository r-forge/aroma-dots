###########################################################################/**
# @RdocClass AbstractProbeSequenceNormalization
#
# @title "The AbstractProbeSequenceNormalization class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a normalization method that corrects for
#  systematic effects in the probe intensities due to differences in
#  probe sequences.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform2".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Requirements}{
#   This class requires that an @see "aroma.core::AromaCellSequenceFile" is 
#   available for the chip type.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AbstractProbeSequenceNormalization", function(...) {
  extend(ProbeLevelTransform2(...), "AbstractProbeSequenceNormalization");
}, abstract=TRUE)



setMethodS3("getTargetFile", "AbstractProbeSequenceNormalization", function(this, ..., verbose=FALSE) {
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


setMethodS3("getAromaCellSequenceFile", "AbstractProbeSequenceNormalization", function(this, ..., force=FALSE) {
  aps <- this$.aps;

  if (force || is.null(aps)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    aps <- AromaCellSequenceFile$byChipType(chipType, ...);
    this$.aps <- aps;
  }

  aps;
}, protected=TRUE)



setMethodS3("fitOne", "AbstractProbeSequenceNormalization", abstract=TRUE, protected=TRUE);

setMethodS3("predictOne", "AbstractProbeSequenceNormalization", abstract=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
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
setMethodS3("process", "AbstractProbeSequenceNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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
  params <- getParameters(this, verbose=less(verbose, 5));

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Get subset to fit
  subsetToFit <- params$subsetToFit;

  # Get subset to update
  subsetToUpdate <- params$subsetToUpdate;

  # Get shift
  shift <- params$shift;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(ds);
  df <- getFile(ds, 1);
  nbrOfCells <- nbrOfCells(df);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);

  hasExcludedMissingSeqs <- FALSE;
  paramsShort <- NULL;
  muT <- NULL;
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
      # Identify probes that cannot be modelled?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!hasExcludedMissingSeqs) {
        verbose && enter(verbose, "Excluding probes  with missing sequences");

        verbose && enter(verbose, "Identifying probes with missing sequences");
        # Locate AromaCellSequenceFile holding probe sequences
        acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));

        missingSeqs <- isMissing(acs, verbose=less(verbose, 5));
        missingSeqs <- whichVector(missingSeqs);
        verbose && cat(verbose, "Cells with unknown sequences:");
        verbose && str(verbose, missingSeqs);
        rm(acs);
        verbose && exit(verbose);

        # Update subset of cell indices for fitting and updating
        subsetToFit <- setdiff(subsetToFit, missingSeqs);
        subsetToUpdate <- setdiff(subsetToUpdate, missingSeqs);
        rm(missingSeqs);

        verbose && cat(verbose, "Cell indices used for fitting:");
        verbose && str(verbose, subsetToFit);
        verbose && cat(verbose, "Cell indices to be updated:");
        verbose && str(verbose, subsetToUpdate);

        # Garbage collection
        gc <- gc();
        verbose && print(verbose, gc);

        hasExcludedMissingSeqs <- TRUE;
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
        verbose && enter(verbose, "Modelling effects for target array");

        verbose && enter(verbose, "Estimating base-count effects for target");
        dfT <- getTargetFile(this, verbose=less(verbose, 5));
        fitT <- fitOne(this, df=dfT, cells=subsetToFit, verbose=less(verbose, 5));
        rm(dfT);
        verbose && print(verbose, fitT);
        verbose && exit(verbose);

        verbose && enter(verbose, "Predicting probe affinities");
        muT <- predictOne(this, fit=fitT, cells=subsetToUpdate, verbose=less(verbose, 5));
        rm(fitT);
        verbose && cat(verbose, "muT:");
        verbose && str(verbose, muT);
        verbose && exit(verbose);

        # Garbage collection
        gc <- gc();
        verbose && print(verbose, gc);

        verbose && exit(verbose);
      } # if (is.null(muT))



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Phase I: Fit base-count effect for the current array
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Getting signals used to fit the model");
      fit <- fitOne(this, df=df, cells=subsetToFit, verbose=less(verbose, 5));
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
      verbose && enter(verbose, "Reading probe signals");
      y <- extractMatrix(df, cells=subsetToUpdate, drop=TRUE);

      # Shift signals?
      if (shift != 0) {
        y <- y + shift;
        verbose && cat(verbose, "Shifted probe signals: ", shift);
      }

      verbose && str(verbose, y);
      verbose && summary(verbose, y);
      verbose && exit(verbose);

      verbose && enter(verbose, "Predicting mean log2 probe signals");
      mu <- predictOne(this, fit=fit, cells=subsetToUpdate, verbose=less(verbose, 5));
      rm(fit);
      verbose && cat(verbose, "mu:");
      verbose && str(verbose, mu);

      verbose && exit(verbose);


      verbose && enter(verbose, "Discrepancy scale factors");
      rho <- (muT-mu);
      rm(mu);
      summary(verbose, rho);
      rho <- 2^rho;
      summary(verbose, rho);

      # Update only subset with "finite" corrections
      keep <- whichVector(is.finite(rho));
      rho <- rho[keep];
      y <- y[keep];
      subsetToUpdateKK <- subsetToUpdate[keep];
      rm(keep);
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Normalizing probe signals");
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
# 2008-07-21
# o The process() is rather generic. Subclasses have to implement fitOne()
#   and predictOne().
# o Created from BaseCountNormalization.R.
############################################################################
