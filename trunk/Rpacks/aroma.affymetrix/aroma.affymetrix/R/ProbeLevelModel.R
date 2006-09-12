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
#   \item{dataSet}{An @see "AffymetrixCelSet" object.}
#   \item{path}{The @character string specifying the path to the directory
#      to contain the parameter-estimate files.}
#   \item{model}{A @character string specifying how PM and MM values
#      should be modelled.  By default only PM signals are used.}
#   \item{...}{Arguments passed to @see "UnitGroupsModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# 
# @author
#
# \seealso{
#   For more details on probe-level models, please see 
#   the \pkg{affyPLM} package.
# }
#*/###########################################################################
setConstructorS3("ProbeLevelModel", function(..., name="modelPLM", model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(UnitGroupsModel(..., name=name), "ProbeLevelModel",
    "cached:.paFile" = NULL,
    "cached:.chipFiles" = NULL,
    "cached:.lastPlotData" = NULL,
    model = model
  )
}, abstract=TRUE)



###########################################################################/**
# @RdocMethod getProbeAffinityClass
#
# @title "Static method to get the ProbeAffinityFile Class object"
#
# \description{
#  @get "title".
#  Any subclass model must provide this method, which should return
#  a @see "R.oo::Class" that subsets the @see "ProbeAffinityFile" class.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "Class" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeAffinityClass", "ProbeLevelModel", abstract=TRUE, static=TRUE);

setMethodS3("getChipEffectClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectFile;
}, static=TRUE, protected=TRUE)

setMethodS3("getChipEffectsClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectSet;
}, static=TRUE, protected=TRUE)



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
setMethodS3("getProbeAffinities", "ProbeLevelModel", function(this, ...) {
  res <- this$.paFile;
  if (!is.null(res))
    return(res);

  # Get the probe-affinity class object
  clazz <- getProbeAffinityClass(this);

  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create probe-affinity file. The CEL set is empty.");

  # Create probe-affinity file from CEL file template
  df <- as.list(ds)[[1]];
  res <- createFrom(df, filename="probeAffinities.CEL", path=getPath(this), ...);

  # Make it into an object of the correct class
  res <- newInstance(clazz, getPathname(res), cdf=getCdf(ds), model=this$model);
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


setMethodS3("readUnits", "ProbeLevelModel", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  # Get what set of probes to read
  stratifyBy <- switch(this$model, pm="pm");

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(cdf, units=units, ..., stratifyBy=stratifyBy);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);

  # Get the CEL intensities by units
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");
  y <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && exit(verbose);
  verbose && str(verbose, y[1]);

  y;
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
#   => For each array and each unit (with G_j groups), we store:
#     G_j theta, G_j sd(theta), G_j isOutlier(theta), i.e. G_j*(float, float, bit)
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

  # Get (and create if missing) the probe-affinity file (one per dataset)
  paf <- getProbeAffinities(this);

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffects(this);

  # Use the same CDF restructor function as the one for the dataset.
  # [This should actually not be required]
  setRestructor(getCdf(ces), getRestructor(cdf));

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
    units <- findUnitsTodo(ces);
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
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
      units <- findUnitsTodo(ces, units=units);
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
  fitfcn <- getFitFunction(this);

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number units per chunk: ", unitsPerChunk);

  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];
    verbose && cat(verbose, "Chunk size: ", length(uu));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the CEL intensities by units
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");
    y <- readUnits(this, units=units[uu], ...);
    verbose && str(verbose, y[1]);
    verbose && exit(verbose);

  
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting unit-group model");
    fit <- lapply(y, FUN=function(unit) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1);
        fitfcn(y);
      })
    })
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && str(verbose, fit[1]);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store probe affinities
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing probe-affinity estimates");
    updateUnits(paf, units=units[uu], data=fit, verbose=verbose);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store chip effects
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing chip-effect estimates");
    updateUnits(ces, units=units[uu], data=fit, verbose=verbose);
    verbose && exit(verbose);

    fit <- NULL; # Not needed anymore

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    gc();
    verbose && exit(verbose);
  } # while()

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
#   \item{verbose}{A @logical or a @see "Verbose" object.}
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
  # Get the chip-effect class object
  clazz <- getChipEffectsClass(this);

  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  cdf <- getCdf(this);
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect file. The CEL set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from dataset");
  ces <- clazz$fromDataSet(dataset=ds, path=getPath(this), verbose=verbose);
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
# 2006-08-28
# o Added plotMvsPosition().
# 2006-08-26
# o Added new getChipEffects().
# 2006-08-25
# o Created from AffymetrixLiWongModel.  So much is in common with the RMA
#   model so a lot can be reused if done cleverly.
############################################################################
