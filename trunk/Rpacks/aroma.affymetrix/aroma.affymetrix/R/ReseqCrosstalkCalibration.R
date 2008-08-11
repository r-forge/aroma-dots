###########################################################################/**
# @RdocClass ReseqCrosstalkCalibration
#
# @title "The ReseqCrosstalkCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the 
#  probe-level signals such that the signals from the four nucleotides 
#  (A, C, G, T) are orthogonal.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{targetAvg}{The signal(s) that the average of the sum of the
#     probe quartets should have after calibration.}
#   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
#     existing indices) used to calculate average in order to rescale to
#     the target average. If @NULL, all probes are considered.}
#   \item{alpha, q, Q}{}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("ReseqCrosstalkCalibration", function(dataSet=NULL, ..., targetAvg=2200, subsetToAvg=NULL, flavor=c("sfit", "expectile"), alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98) {
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

    # Argument 'targetAvg':
    if (!is.null(targetAvg)) {
      targetAvg <- Arguments$getDouble(targetAvg, range=c(0, Inf));
    }
  
    # Argument 'subsetToAvg':
    if (is.null(subsetToAvg)) {
    } else if (is.character(subsetToAvg)) {
      if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToAvg': ", subsetToAvg);
      }
      extraTags <- c(extraTags, subsetToAvg=subsetToAvg);
    } else {
      subsetToAvg <- Arguments$getIndices(subsetToAvg, 
                                          range=c(1, nbrOfCells(cdf)));
      subsetToAvg <- unique(subsetToAvg);
      subsetToAvg <- sort(subsetToAvg);
    }
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  

  extend(ProbeLevelTransform(dataSet=dataSet, ...), "ReseqCrosstalkCalibration",
    .targetAvg = targetAvg,
    .subsetToAvg = subsetToAvg,
    .extraTags = extraTags,
    .flavor = flavor,
    .alpha = alpha,
    .q = q,
    .Q = Q
  )
})


setMethodS3("clearCache", "ReseqCrosstalkCalibration", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".subsetToAvgExpanded")) {
    this[[ff]] <- NULL;
  }

  # Then for this object 
  NextMethod("clearCache", object=this, ...);
})


setMethodS3("getAsteriskTags", "ReseqCrosstalkCalibration", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)



setMethodS3("getSubsetToAvg", "ReseqCrosstalkCalibration", function(this, ..., verbose=FALSE) {
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

  # Expand?
  if (is.character(subsetToAvg)) {
    if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToAvg: ", subsetToAvg);

      # Look up in cache
      subset <- this$.subsetToAvgExpanded;
      if (is.null(subset)) {
        dataSet <- getInputDataSet(this);
        cdf <- getCdf(dataSet);
  
        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);
  
        # Identify units to be excluded
        if (subsetToAvg == "-X") {
          subset <- getUnitsOnChromosome(gi, 23, .checkArgs=FALSE);
        } else if (subsetToAvg == "-Y") {
          subset <- getUnitsOnChromosome(gi, 24, .checkArgs=FALSE);
        } else if (subsetToAvg == "-XY") {
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
        this$.subsetToAvgExpanded <- subset;
      }

      subsetToAvg <- subset;
      rm(subset);

      verbose && exit(verbose);
    }
  }

  subsetToAvg;
}, protected=TRUE);



setMethodS3("getParameters", "ReseqCrosstalkCalibration", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, expand=expand, ...);

  params <- c(params, list(
    targetAvg = this$.targetAvg,
    subsetToAvg = this$.subsetToAvg,
    flavor = this$.flavor,
    alpha = this$.alpha,
    q = this$.q,
    Q = this$.Q
  ));

  # Expand?
  if (expand) {
    params$subsetToAvg <- getSubsetToAvg(this);
  }

  params;
}, private=TRUE)



setMethodS3("rescaleByAll", "ReseqCrosstalkCalibration", function(this, yAll, params, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  targetAvg <- params$targetAvg;
  nt <- length(targetAvg);
  if (nt != 1) {
    throw("In order rescale towards a global target average ('rescaleBy' == \"all\"), argument 'targetAvg' must be a scalar: ", paste(targetAvg, collapse=","));
  }
 
  if (verbose) {
    enter(verbose, "Rescaling toward target average");
    cat(verbose, "Target average: ", targetAvg);
    if (!is.null(params$subsetToAvg)) {
      cat(verbose, "Using subset of cells for estimate of target average:");
      src <- attr(params$subsetToAvg, "src");
      if (is.null(src))
         src <- params$subsetToAvg;
      str(verbose, src);
    }
    cat(verbose, "yAll: ");
    str(verbose, yAll);
  }

  # Total number of values
  n0 <- length(yAll);

  # Average of *all* values
  yAvg0 <- median(yAll, na.rm=TRUE);

  if (!is.null(params$subsetToAvg)) {
    y <- yAll[params$subsetToAvg];
    n <- length(y);
    if (n == 0) {
      throw("Cannot rescale to target average. There are no cells to average over.");
    }
    yAvg <- median(y, na.rm=TRUE);
    verbose && printf(verbose, "yAvg (using %d/%.1f%% summed pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
    rm(y, n);
  } else {
    yAvg <- yAvg0;
    verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
  }
    
  if (!is.finite(yAvg)) {
    throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);
  }

  b <- targetAvg/yAvg;
  verbose && printf(verbose, "Scale factor: %.2f\n", b);
  yAll <- b*yAll;
  fit <- list(all=list(b=b));
  verbose && exit(verbose);

  attr(yAll, "fit") <- fit;

  verbose && cat(verbose, "Rescaling parameter estimates:");
  verbose && str(verbose, fit);
  verbose && exit(verbose);

  yAll;
}, protected=TRUE)




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
setMethodS3("process", "ReseqCrosstalkCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calibrating data set for cross talk");

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
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # To be retrieved when needed.
  cellQuartets <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "process.ReseqCrosstalkCalibration";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate some model fit parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Compressing model parameter to a short format");
  paramsShort <- params;
  paramsShort$subsetToAvg <- NULL;
#  paramsShort$subsetToAvgIntervals <- seqToIntervals(params$subsetToAvg);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- nbrOfArrays(ds);
  verbose && enter(verbose, "Calibrating ", nbrOfArrays, " arrays");
  verbose && cat(verbose, "Path: ", outputPath);
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
      if (is.null(cellQuartets)) {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Identify the cell-index matrix
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        verbose && enter(verbose, "Identifying cell-index matrix");
        cellQuartets <- getCellQuartets(cdf, verbose=verbose);
#        cellQuartets <- cellQuartets[1:1000,];
        verbose && cat(verbose, "cellQuartets:");
        verbose && str(verbose, cellQuartets);
        dim <- dim(cellQuartets);
        gc <- gc();
        verbose && print(verbose, gc);
        verbose && exit(verbose);
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities");
      yAll <- getData(df, fields="intensities", ...)$intensities;
      verbose && str(verbose, yAll);
      verbose && exit(verbose);
    
      modelFit <- list(
        paramsShort=paramsShort
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fitting crosstalk model
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Fitting crosstalk model");
      verbose && cat(verbose, "Number of dimensions: ", dim[2]);
      verbose && cat(verbose, "Number of data points: ", dim[1]);
  
      verbose && enter(verbose, "Extracting data for fitting model");
      y <- yAll[cellQuartets];
      dim(y) <- dim;
      verbose && str(verbose, y);
      verbose && exit(verbose);

      verboseL <- (verbose && isVisible(verbose, -50));

      verbose && enter(verbose, "Fitting");
      verbose && cat(verbose, "Flavor: ", params$flavor);
      fit <- fitMultiDimensionCone(y, flavor=params$flavor,
                          alpha=params$alpha, q=params$q, Q=params$Q, 
                                                     verbose=verboseL);
      verbose && print(verbose, fit, level=-5);
      verbose && str(verbose, fit);
      verbose && exit(verbose);

      callHooks(sprintf("%s.onFitOne", hookName), df=df, y=y, fit=fit, ...);
      rm(y); # Not needed anymore

      verbose && exit(verbose);

      # Store crosstalk model fit
      modelFit$rccFit <- fit;


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Backtransforming (calibrating)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Calibrating data");

      verbose && enter(verbose, "Extracting data to be calibrated");
      y <- yAll[cellQuartets];
      dim(y) <- dim;
      verbose && str(verbose, y);
      verbose && exit(verbose);

      verbose && enter(verbose, "Backtransforming");
      yC <- backtransformMultiDimensionCone(y, fit=fit);
      rm(y, fit); # Not needed anymore
      verbose && str(verbose, yC);
      yAll[cellQuartets] <- yC;
      rm(yC); # Not needed anymore
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Rescaling toward target average?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!is.null(params$targetAvg)) {
        yAll <- rescaleByAll(this, yAll=yAll, params=params, 
                           setsOfProbes=setsOfProbes, verbose=less(verbose));
        fit <- attr(yAll, "fit");
        fit$params <- NULL;
        fit$paramsShort <- paramsShort;
        modelFit$rescaleFit <- fit;
        rm(fit);

        # Garbage collect
        gc <- gc();
        verbose && print(verbose, gc);
      }


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store model fit 
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname);
      fitPathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      saveObject(modelFit, file=fitPathname);
      verbose && str(verbose, modelFit, level=-50);
      rm(modelFit);


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

    # Test if calibrated data file can be retrieved
    dfC <- newInstance(df, pathname);

    rm(df); # Not needed anymore

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  # Garbage collect
  rm(ds, cellQuartets);

  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
}) # process()





############################################################################
# HISTORY:
# 2008-08-10
# o Created from AllelicCrosstalkCalibration.R.
############################################################################
