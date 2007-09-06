###########################################################################/**
# @RdocClass AllelicCrosstalkCalibration
#
# @title "The AllelicCrosstalkCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the 
#  probe-level signals such that the signals from the two alleles are 
#  orthogonal.
#  The method fits and calibrates PM signals only.  MM signals will not 
#  affect the model fitting and are unaffected.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{targetAvg}{The signal(s) that either the average of the sum
#     (if one target value) or the average of each of the alleles
#     (if two target values) should have after calibration.}
#   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
#     existing indices) used to calculate average in order to rescale to
#     the target average. If @NULL, all probes are considered.}
#   \item{alpha, q, Q}{}
# }
#
# \section{What probe signals are updated?}{ 
#   Calibration for crosstalk between allele signals applies by definition
#   only SNP units.  It is only PM probes that will be calibrated.
#   Note that, non-calibrated signals will be saved in the output files.
# }
#
# \section{What probe signals are used to fit model?}{
#   All PM probe pairs are used to fit the crosstalk model.
#   In the second step where signals are rescaled to a target average,
#   it is possible to specify the set of cells that should be included
#   when estimating the target average.
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AllelicCrosstalkCalibration", function(..., targetAvg=c(2200, 2200), subsetToAvg=NULL, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98) {
  if (!is.null(targetAvg)) {
    targetAvg <- Arguments$getDoubles(targetAvg, range=c(0, Inf));
    if (!length(targetAvg) %in% 1:2) {
      throw("Argument 'targetAvg' must be of length one or two: ", length(targetAvg));
    }
  }

  extend(ProbeLevelTransform(...), "AllelicCrosstalkCalibration",
    .targetAvg = targetAvg,
    .subsetToAvg = subsetToAvg,
    .alpha = alpha,
    .q = q,
    .Q = Q
  )
})


setMethodS3("getSubsetToAvg", "AllelicCrosstalkCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  subsetToAvg <- this$.subsetToAvg;
  if (FALSE) {
    this$.subsetToAvg <- subsetToAvg;
  }

  subsetToAvg;
}, protected=TRUE);



setMethodS3("getParameters", "AllelicCrosstalkCalibration", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  params <- c(params, list(
    targetAvg = this$.targetAvg,
    subsetToAvg = getSubsetToAvg(this),
    alpha = this$.alpha,
    q = this$.q,
    Q = this$.Q
  ));

  params;
}, private=TRUE)


setMethodS3("rescale", "AllelicCrosstalkCalibration", function(this, yAll, params, setsOfProbes, ..., verbose=FALSE) {
  nbrOfPairs <- length(setsOfProbes);
  nt <- length(params$targetAvg);
  if (nt == 1) {
    method <- "sum";
  } else if (nt == 2) {
    method <- "allele";
  }

  if (verbose) {
    enter(verbose, "Rescaling toward target average");
    cat(verbose, "Target average(s): %s", paste(params$targetAvg, collapse=", "));
    printf(verbose, "Method: %s\n", method);
    if (!is.null(params$subsetToAvg)) {
      cat(verbose, "Using subset of cells for estimate of target average:");
      str(verbose, params$subsetToAvg);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale based on y = yA+yB
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "sum") {
    for (kk in seq_len(nbrOfPairs)) {
      name <- names(setsOfProbes)[kk];
      verbose && enter(verbose, sprintf("Allele basepair #%d ('%s') of %d", kk, name, nbrOfPairs));
  
      # Get data pairs
      idxAB <- setsOfProbes[[name]];
      idxAB <- matrix(idxAB, ncol=2);

      # Sum y=yA+yB
      y <- yAll[idxAB[,1]]+yAll[idxAB[,2]];
      n <- length(y);
      n0 <- n;
    
      # Calculate current average
      yAvg <- median(y, na.rm=TRUE);
      yAvg0 <- yAvg;
      rm(y);
    
      if (!is.null(params$subsetToAvg)) {
        keep <- matrix((idxAB %in% params$subsetToAvg), ncol=2);
        keep <- (keep[,1] & keep[,2]);
        idxAB <- idxAB[keep,,drop=FALSE];
        rm(keep);

        # Sum y=yA+yB
        y <- yAll[idxAB[,1]]+yAll[idxAB[,2]];
        n <- length(y);
    
        if (n == 0) {
          throw("Cannot rescale to target average. After taking the intersect of the subset of cells to be used, there are no cells left.");
        }
    
        yAvg <- median(y, na.rm=TRUE);
        rm(y);
    
        verbose && printf(verbose, "yAvg (using %d/%.1f%% summed pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
      } else {
        verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
      }
    
      if (!is.finite(yAvg))
        throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);
    
      # Rescale
      b <- params$targetAvg/yAvg;
      verbose && printf(verbose, "scale factor: %.2f\n", b);
    
      idxAB <- setsOfProbes[[name]];
      yAll[idxAB] <- b*yAll[idxAB];
  
      rm(idx);
      verbose && exit(verbose);
    } # for (kk in ...)
    verbose && exit(verbose);
  } # if (method == "sum")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale based on (yA,yB)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "allele") {
    for (kk in seq_len(nbrOfPairs)) {
      name <- names(setsOfProbes)[kk];
      verbose && enter(verbose, sprintf("Allele basepair #%d ('%s') of %d", kk, name, nbrOfPairs));
  
      # Get data pairs
      idxAB <- setsOfProbes[[name]];
      idxAB <- matrix(idxAB, ncol=2);

      # For each allele
      for (cc in 1:2) {
        idx <- idxAB[,cc];
        y <- yAll[idx];
        n <- length(y);
        n0 <- n;
    
        # Calculate current average
        yAvg <- median(y, na.rm=TRUE);
        yAvg0 <- yAvg;
        rm(y);
    
        if (!is.null(params$subsetToAvg)) {
          idx <- idxAB[,cc];
          idx <- intersect(idx, params$subsetToAvg);
          y <- yAll[idx];
          n <- length(y);
    
          if (n == 0) {
            throw("Cannot rescale to target average. After taking the intersect of the subset of cells to be used, there are no cells left.");
          }
    
          yAvg <- median(y, na.rm=TRUE);
          rm(y);
    
          verbose && printf(verbose, "yAvg (using %d/%.1f%% pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
        } else {
          verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
        }
    
        if (!is.finite(yAvg))
          throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);
    
        # Rescale
        b <- params$targetAvg[cc]/yAvg;
        verbose && printf(verbose, "scale factor: %.2f\n", b);
    
        idx <- idxAB[,cc];
        yAll[idx] <- b*yAll[idx];
  
        rm(idx);
      } # for (cc ...)
      verbose && exit(verbose);
    } # for (kk in ...)
    verbose && exit(verbose);
  } # if (method == "allelic")

  yAll;
}, protected=TRUE)


setMethodS3("getDataPairs", "AllelicCrosstalkCalibration", function(this, array, cs=NULL, ..., verbose=FALSE) {
  if (is.null(cs)) {
    cs <- getInputDataSet(this);
  }
  cdf <- getCdf(cs);

  verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
  setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading all probe intensities");
  cf <- getFile(cs, array);
  yAll <- getData(cf, fields="intensities", ...)$intensities;
  verbose && exit(verbose);

  nbrOfPairs <- length(setsOfProbes);
  res <- vector("list", nbrOfPairs);
  names(res) <- names(setsOfProbes);

  verbose && enter(verbose, "Extracting data pairs");
  for (kk in seq_len(nbrOfPairs)) {
    name <- names(setsOfProbes)[kk];
    basepair <- unlist(strsplit(name, split=""));
    idx <- setsOfProbes[[name]];
    y <- matrix(yAll[idx], ncol=2);
    colnames(y) <- c("A", "B");
    res[[kk]] <- y;
    rm(y);
  }
  verbose && exit(verbose);

  res;
}, protected=TRUE);


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
setMethodS3("process", "AllelicCrosstalkCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calibrates data set for allelic cross talk");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already calibrated");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require("sfit") || throw("Package not loaded: sfit.");

  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # To be retrieved when needed.
  setsOfProbes <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "process.AllelicCrosstalkCalibration";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- nbrOfArrays(ds);
  verbose && enter(verbose, "Calibrating ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);
  dataFiles <- list();
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already calibrated?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Calibrated data file already exists: ", pathname);
    } else {
      if (is.null(setsOfProbes)) {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Identify the cell indices for each possible allele basepair.
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
        setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
        gc <- gc();
        verbose && print(verbose, gc);
        verbose && exit(verbose);
      }
  
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities");
      yAll <- getData(df, fields="intensities", ...)$intensities;
      verbose && exit(verbose);
    
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Calibrating
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      callHooks(sprintf("%s.onBegin", hookName), df=df, setsOfProbes=setsOfProbes, ...);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fitting
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nbrOfPairs <- length(setsOfProbes);
      fits <- vector("list", nbrOfPairs);
      names(fits) <- names(setsOfProbes);
      verbose && enter(verbose, "Fitting calibration model");
      for (kk in seq_len(nbrOfPairs)) {
        name <- names(setsOfProbes)[kk];
        verbose && enter(verbose, sprintf("Allele basepair #%d ('%s') of %d", kk, name, nbrOfPairs));
        basepair <- unlist(strsplit(name, split=""));
        idx <- setsOfProbes[[name]];
    
        verbose && enter(verbose, "Fitting");
        y <- matrix(yAll[idx], ncol=2);
#        callHooks(sprintf("%s.onData", hookName), df=df, y=y, basepair=basepair, ...);
        fits[[kk]] <- fitGenotypeCone(y, alpha=params$alpha, q=params$q, Q=params$Q);
        verbose && print(verbose, fits[[kk]], level=-5);
        verbose && exit(verbose);

        rm(y, idx); # Not needed anymore
        gc <- gc();

#        callHooks(sprintf("%s.onFit", hookName), df=df, fit=fits[[kk]], ...);
        verbose && exit(verbose);
      } # for (kk in ...)
      verbose && exit(verbose);


      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);
      modelFit <- list(
        params=getParameters(this),
        fits=fits
      );
      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      rm(modelFit);

      gc <- gc();
      verbose && print(verbose, gc);
    

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Harmonizing parameter estimates?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Here we can harmonize the estimates, e.g. make all offset estimates
      # the same.  This might be useful if we want to correct other probes
      # not included above such as CN probes on SNP 6.0. /HB 2007-09-05


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Backtransforming (calibrating)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Backtransforming (calibrating) data");
      for (kk in seq_len(nbrOfPairs)) {
        name <- names(setsOfProbes)[kk];
        verbose && enter(verbose, sprintf("Allele basepair #%d ('%s') of %d", kk, name, nbrOfPairs));

        idx <- setsOfProbes[[name]];
        y <- matrix(yAll[idx], ncol=2);
        yC <- backtransformGenotypeCone(y, fit=fits[[kk]]);
        yAll[idx] <- yC;
    
#        callHooks(sprintf("%s.onUpdated", hookName), df=df, y=y, basepair=basepair, fit=fits[[kk]], yC=yC,...);
        rm(idx, y, yC);
        gc <- gc();

        verbose && exit(verbose);
      } # for (kk in ...)
      verbose && exit(verbose);

      # Not needed anymore
      rm(fits);
      gc <- gc();

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Rescaling toward target average?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!is.null(params$targetAvg)) {
        yAll <- rescale(this, yAll=yAll, params=params, setsOfProbes=setsOfProbes, verbose=less(verbose));
      }



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing calibrated data");
    
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

    # Retrieving calibrated data file
    dfC <- newInstance(df, pathname);

    # CDF inheritance
    setCdf(dfC, cdf);

#    callHooks(sprintf("%s.onExit", hookName), df=df, dfC=dfC, ...);

    # Record
    dataFiles[[kk]] <- dfC;

    rm(df, dfC);

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  # Garbage collect
  rm(dataFiles, ds, setsOfProbes);
  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})




############################################################################
# HISTORY:
# 2007-09-05
# o Now the rescaling can be done either on (yA,yB) separately or on 
#   y=yA+yB.  If targetAvg has two values the former, otherwise the latter.
# o Now AllelicCrosstalkCalibration recognizes argument 'subsetToAvg'.
# o Now process() stores the crosstalk settings and estimated parameters
#   to file. May be useful if one wants to go back and look at the details.
#   One day we might get around to store this information in the CEL file
#   headers.
# o Now process() first fits the crosstalk model for all basepairs, then
#   backtransform the signals, then optional rescale signals to target
#   average, then saves the calibrated signals.
# o SPEED UP: Now getAlleleProbePairs() is only called if data needs to be
#   calibrated, i.e. if already calibrated it is not loaded.
# o CLEAN UP: Now the code of calibrateAllelicCrosstalk() is included here.
# 2007-03-29
# o Now 'targetAvg' defaults to 2200 so that allele A and allele B signals
#   are rescaled to be one the same scale.  If so,  it does not make sense
#   to do background correction afterwards.
# o Added getParameters().
# o Added support for arguments 'targetAvg', 'alpha', 'q', and 'Q'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessing class.
# o Now this pre-processor output results to probeData/.
# o Renamed from AllelicCrosstalkCalibrator.
# 2006-11-18
# o Removed version and subversion tags, and related functions. 
#   Now getTags() returns the tags of the input data set plus any tags 
#   of this instance.
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
