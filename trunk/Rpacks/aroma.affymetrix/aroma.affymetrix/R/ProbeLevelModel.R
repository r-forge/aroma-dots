###########################################################################/**
# @RdocClass ProbeLevelModel
#
# @title "The ProbeLevelModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents any probe-level model (PLM).
#  To quote the \pkg{affyPLM} package:
#    "A [...] PLM is a model that is fit to probe-intensity data. 
#     More specifically, it is where we fit a model with probe level
#     and chip level parameters on a probeset by probeset basis",
#  where the more general case for a probeset is a unit group
#  using Affymetrix CDF terms.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "UnitGroupsModel".}
#   \item{probeModel}{A @character string specifying how PM and MM values
#      should be modelled.  By default only PM signals are used.}
#   \item{standardize}{If @TRUE, chip-effect and probe-affinity estimates are
#      rescaled such that the product of the probe affinities is one.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   In order to minimize the risk for mistakes, but also to be able compare
#   results from different PLMs, all PLM subclasses must return parameter
#   estimates that meet the following requirements:
#   \enumerate{
#     \item All parameter estimates must be (stored and returned) on the
#       intensity scale, e.g. log-additive models such as @see "RmaPlm"
#       have to transform the parameters on the log-scale to the intensity
#       scale.
#     \item The probe-affinity estimates \eqn{\phi_j} for a unit group
#       must meet the constraint such that \eqn{\prod_j \phi_j = 1},
#       or equivalently \eqn{\sum_j \log(\phi_j) = 0}.
#   }
#   Note that the above probe-affinity constraint guarantees that the
#   estimated chip effects across models are on the same scale.
# }
# 
# @author
#
# \seealso{
#   For more details on probe-level models, please see 
#   the \pkg{affyPLM} package.
# }
#*/###########################################################################
setConstructorS3("ProbeLevelModel", function(..., tags=NULL, probeModel=c("pm", "mm", "pmmm"), standardize=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  if (probeModel != "pm") {
    tags <- c(toupper(probeModel), tags);
  }

  extend(UnitGroupsModel(..., tags=tags), "ProbeLevelModel",
    "cached:.paFile" = NULL,
    "cached:.chipFiles" = NULL,
    "cached:.lastPlotData" = NULL,
    probeModel = probeModel,
    standardize=standardize
  )
}, abstract=TRUE)


setMethodS3("getChipEffectSetClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectSet;
}, static=TRUE);


setMethodS3("getRootPath", "ProbeLevelModel", function(this, ...) {
  "plmData";
})


###########################################################################/**
# @RdocMethod getProbeAffinities
#
# @title "Gets the probe affinities for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "ProbeAffinityFile" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeAffinities", "ProbeLevelModel", function(this, ..., .class=ProbeAffinityFile) {
  res <- this$.paFile;
  if (!is.null(res))
    return(res);

  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create probe-affinity file. The CEL set is empty.");

  # Create probe-affinity file from CEL file template
  df <- as.list(ds)[[1]];
  res <- createFrom(df, filename="probeAffinities.cel", path=getPath(this), ...);

  # Make it into an object of the correct class
  clazz <- .class;
  res <- newInstance(clazz, getPathname(res), cdf=getCdf(ds), probeModel=this$probeModel);
  this$.paFile <- res;

  res;
})



###########################################################################/**
# @RdocMethod getFitFunction
#
# @title "Static method to get the low-level function that fits the PLM"
#
# \description{
#  @get "title".
#  Any subclass model must provide this method, which should return
#  a @function that accepts an IxJ @matrix.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @function.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "ProbeLevelModel", abstract=TRUE, static=TRUE);


setMethodS3("getFitUnitFunction", "ProbeLevelModel", function(this, ...) {
  # Get the fit function for a single set of intensities
  fitfcn <- getFitFunction(this, ...);

  # Create the one for all blocks in a unit
  fitUnit <- function(unit, ...) {
    lapply(unit, FUN=function(group) {
      y <- .subset2(group, 1);
      fitfcn(y);
    })
  }

  fitUnit;
})





###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset) 
#  specially structured for this PLM.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Arguments passed to \code{getCellIndices()} of the 
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the \code{readUnits()} method of the 
#     @see "AffymetrixCelSet" class.}
# }
#
# \value{
#  Returns the @list structure that \code{readUnits()} of 
#  @see "AffymetrixCelSet" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "ProbeLevelModel", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ...);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);

  # Get the CEL intensities by units
  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");
  res <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && exit(verbose);
  verbose && str(verbose, res[1]);

  res;
})


setMethodS3("getCellIndices", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get what set of probes to read
  stratifyBy <- switch(this$probeModel, mm="mm", pm="pm", pmmm="pmmm");

  # Get the CDF cell indices
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  verbose && enter(verbose, "Identifying CDF cell indices");
<<<<<<< .mine
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);
  verbose && str(verbose, list(...));
=======
  verbose && cat(verbose, "Stratify by: ", stratifyBy);
>>>>>>> .r1300
  cells <- getCellIndices(cdf, ..., stratifyBy=stratifyBy);
  verbose && exit(verbose);
  
  cells;
})


###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies units for which the PLM has still not been fitted to"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   Internally this methods calls the same method for the 
#   @see "ChipEffectSet" class.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ces <- getChipEffects(this, verbose=verbose);
  findUnitsTodo(ces, ..., verbose=verbose);
})



###########################################################################/**
# @RdocMethod fit
#
# @title "Fits a unit groups (probeset) model"
#
# \description{
#  @get "title" to all or to a subset of the units.
#  All estimates are stored to file.
#  The non-array specific parameter estimates together with standard deviation
#  estimates and convergence information are stored in one file.
#  The parameter estimates specific to each array, typically "chip effects", 
#  are stored in array specific files.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the indices of the units fitted, or @NULL if no units had to be fitted.
# }
#
# \details{
#   Dataset-specific estimates [L = nbr of probes]:
#    phi [L doubles] (probe affinities), sd(phi) [L doubles], 
#    isOutlier(phi) [L logicals]
#
#   Algorithm-specific results:
#    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
#    dTheta [1 double]
#    sd(eps) - [1 double] estimated standard deviation of the error term
#
#   Array-specific estimates [K = nbr of arrays]:
#    theta [K doubles] (chip effects), sd(theta) [K doubles], 
#    isOutlier(theta) [K logicals]
#   
#   For each array and each unit group, we store:
#     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
#   => For each array and each unit (with \eqn{G_j} groups), we store:
#     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta), 
#   i.e. \eqn{G_j}*(float, float, bit).
#   For optimal access we store all thetas first, then all sd(theta) and the
#   all isOutlier(theta).
#   To keep track of the number of groups in each unit, we have to have a
#   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
#   AffymetrixCdfFile class.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "ProbeLevelModel", function(this, units="remaining", ..., unitsPerChunk=moreUnits*100000/length(getDataSet(this)), moreUnits=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting model of class ", class(this)[1], ":");

  verbose && print(verbose, this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, verbose=less(verbose));
    nbrOfUnits <- length(units);
    verbose && exit(verbose);
  } else {
    # Fit only unique units
    units <- unique(units);
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }

  # Nothing more to do?
  if (nbrOfUnits == 0)
    return(NULL);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit the model in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  # Get (and create if missing) the probe-affinity file (one per dataset)
  paf <- getProbeAffinities(this, verbose=less(verbose));

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffects(this, verbose=less(verbose));

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number units per chunk: ", unitsPerChunk);

  # Time the fitting.
  startTime <- processTime();

  timers <- list(total=0, read=0, fit=0, writePaf=0, writeCes=0, gc=0);

  count <- 1;
  while (length(idxs) > 0) {
    tTotal <- processTime();

    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the CEL intensities by units
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nbrOfArrays <- length(ds);
    verbose && enter(verbose, "Reading probe intensities from ", nbrOfArrays, " arrays");
    tRead <- processTime();
    y <- readUnits(this, units=units[uu], ..., verbose=less(verbose));
    timers$read <- timers$read + (processTime() - tRead);
    verbose && str(verbose, y[1]);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting unit-group model");
    tFit <- processTime();
    fit <- lapply(y, FUN=fitUnit);
    timers$fit <- timers$fit + (processTime() - tFit);
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && str(verbose, fit[1]);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store probe affinities
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing probe-affinity estimates");
    tWritePaf <- processTime();
    updateUnits(paf, units=units[uu], data=fit, verbose=less(verbose));
    timers$writePaf <- timers$writePaf + (processTime() - tWritePaf);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store chip effects
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing chip-effect estimates");
    tWriteCes <- processTime();
    updateUnits(ces, units=units[uu], data=fit, verbose=less(verbose));
    timers$writeCes <- timers$writeCes + (processTime() - tWriteCes);
    verbose && exit(verbose);

    fit <- NULL; # Not needed anymore

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc();
    timers$gc <- timers$gc + (processTime() - tGc);

    timers$total <- timers$total + (processTime() - tTotal);

    verbose && exit(verbose);
  } # while()

  totalTime <- processTime() - startTime;
  if (verbose) {
    nunits <- length(units);
    t <- totalTime[3];
    printf(verbose, "Total time for all units across all %d arrays: %.2fs == %.2fmin\n", nbrOfArrays, t, t/60);
    t <- totalTime[3]/nunits
    printf(verbose, "Total time per unit across all %d arrays: %.2fs/unit\n", nbrOfArrays, t);
    t <- totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time per unit and array: %.3gms/unit & array\n", 1000*t);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", nbrOfUnits(cdf), t/60, t/3600);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits;
    printf(verbose, "Total time for complete dataset: %.2fmin = %.2fh\n", t/60, t/3600);
    # Get distribution of what is spend where
    timers$write <- timers$writePaf + timers$writeCes;
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for writing chip-effects), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeCes"]/t["write"], t["gc"]);
  }

  invisible(units);
})




###########################################################################/**
# @RdocMethod getChipEffects
#
# @title "Gets the chip-effect files for this model"
#
# \description{
#  @get "title".  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#  Returns a @see "ChipEffectSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffects", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ces <- this$.ces;
  if (!is.null(ces))
    return(ces);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect file. The CEL set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from dataset");
  # Gets the ChipEffects Class object
  clazz <- getChipEffectSetClass(this);
  ces <- clazz$fromDataSet(dataset=ds, path=getPath(this), verbose=less(verbose));
  verbose && exit(verbose);

  # Store in cache
  this$.ces <- ces;

  ces;
})



setMethodS3("plotMvsPosition", "ProbeLevelModel", function(this, sample, ..., annotate=TRUE) {
  ces <- getChipEffects(this);
  ce <- getFile(ces, sample);
  cesAvg <- getAverageFile(ces);
  res <- plotMvsPosition(ce, reference=cesAvg, ..., annotate=annotate);  
  if (annotate) {
    stext(getLabel(this), side=1, pos=1, line=-1, cex=0.7, col="darkgrey");
  }

  this$lastPlotData <- res;
  this$lastPlotSample <- sample;

  invisible(res);
})


setMethodS3("highlight", "ProbeLevelModel", function(this, ...) {
  ces <- getChipEffects(this);
  ce <- getFile(ces, this$lastPlotSample);
  highlight(ce, ...);
})



############################################################################
# HISTORY:
# 2006-09-26
# o Fixed the timers for fit(). They only worked so and so before and only
#   only Windows.  Using processTime()
# 2006-09-14
# o Added detailed timing information to the verbose output of fit().
# 2006-09-11
# o Added argument ProbeLevelModel(..., standardize=TRUE) to make the 
#   results from different PLMs be on the same scale.
# 2006-09-10
# o Added findUnitsTodo().
# o Update fit() to make use of new ChipEffects and ProbeAffinity classes.
#   The code is much cleaner now!
# 2006-08-28
# o Added plotMvsPosition().
# 2006-08-26
# o Added new getChipEffects().
# 2006-08-25
# o Created from AffymetrixLiWongModel.  So much is in common with the RMA
#   model so a lot can be reused if done cleverly.
############################################################################
