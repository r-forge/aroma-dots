###########################################################################/**
# @RdocClass ProbeLevelModel
#
# @title "The ProbeLevelModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a probe-level model (PLM) as defined
#  by the \pkg{affyPLM} package:
#    "A [...] PLM is a model that is fit to probe-intensity data. 
#     More specifically, it is where we fit a model with probe level
#     and chip level parameters on a probeset by probeset basis",
#  where the more general case for a probeset is a \emph{unit group}
#  in Affymetrix CDF terms.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "UnitModel".}
#   \item{tags}{A @character @vector of tags to be added.}
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
#   results from different PLMs, all PLM subclasses must meet the following
#   criteria:
#   \enumerate{
#     \item All parameter estimates must be (stored and returned) on the
#       intensity scale, e.g. log-additive models such as @see "RmaPlm"
#       have to transform the parameters on the log-scale to the intensity
#       scale.
#     \item The probe-affinity estimates \eqn{\phi_j} for a unit group
#       must be constrained such that \eqn{\prod_j \phi_j = 1},
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
setConstructorS3("ProbeLevelModel", function(..., tags=NULL, probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm"), standardize=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  
  # Add tags
  if (probeModel != "pm") {
    tags <- c(toupper(probeModel), tags);
  }

  extend(UnitModel(..., tags=tags, parSet=list(probeModel=probeModel)), "ProbeLevelModel",
    "cached:.paf" = NULL,
    "cached:.ces" = NULL,
    "cached:.rs" = NULL,
    "cached:.ws" = NULL,
    "cached:.lastPlotData" = NULL,
    probeModel = probeModel,
    standardize = standardize,
    shift = 0
  )
}, abstract=TRUE)

setMethodS3("clearCache", "ProbeLevelModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".paf", ".ces", ".rs", ".lastPlotData")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("getParameterSet", "ProbeLevelModel", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$probeModel <- this$probeModel;
  params$shift <- this$shift;
  params;
}, private=TRUE)



setMethodS3("getRootPath", "ProbeLevelModel", function(this, ...) {
  "plmData";
}, private=TRUE)


###########################################################################/**
# @RdocMethod getProbeAffinityFile
# @aliasmethod getProbeAffinities
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
#   \item{.class}{A @see "ProbeAffinityFile" \emph{class}.}
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
setMethodS3("getProbeAffinityFile", "ProbeLevelModel", function(this, ..., .class=ProbeAffinityFile) {
  paf <- this$.paf;
  if (!is.null(paf))
    return(paf);

  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create probe-affinity file. There are no CEL files in the data set.");

  # Create probe-affinity file from CEL file template
  df <- getFile(ds, 1);
  paf <- createFrom(df, filename="probeAffinities.cel", path=getPath(this), 
                                                         method="create", ...);

  # Make it into an object of the correct class
  paf <- newInstance(.class, getPathname(paf), cdf=getCdf(ds), 
                                                  probeModel=this$probeModel);

  this$.paf <- paf;

  paf;
})



setMethodS3("getProbeAffinities", "ProbeLevelModel", function(this, ...) {
  getProbeAffinityFile(this, ...);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getChipEffectSet
# @aliasmethod getChipEffects
#
# @title "Gets the set of chip effects for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
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
setMethodS3("getChipEffectSet", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
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
    throw("Cannot create chip-effect set. The CEL set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from data set");
  # Gets the ChipEffects Class object
  clazz <- getChipEffectSetClass(this);
  ces <- clazz$fromDataSet(dataSet=ds, path=getPath(this), 
                                                    verbose=less(verbose));
  verbose && exit(verbose);

  # Store in cache
  this$.ces <- ces;

  ces;
})

setMethodS3("getChipEffects", "ProbeLevelModel", function(this, ...) {
  getChipEffectSet(this, ...);
}, protected=TRUE)


setMethodS3("getChipEffectSetClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectSet;
}, static=TRUE, private=TRUE)



setMethodS3("getResidualSet", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  rs <- this$.rs;
  if (!is.null(rs))
    return(rs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create residuals files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create residuals set. The data set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from data set");
  # Gets the ResidualSet Class object
  clazz <- getResidualSetClass(this);
  rs <- clazz$fromDataSet(dataSet=ds, path=getPath(this), 
                                                     verbose=less(verbose));
  # make sure CDF is inherited
  setCdf(rs, getCdf(ds));
  verbose && exit(verbose);

  # Store in cache
  this$.rs <- rs;

  rs;
})

setMethodS3("getResidualSetClass", "ProbeLevelModel", function(static, ...) {
  ResidualSet;
}, static=TRUE, private=TRUE)


setMethodS3("getWeightsSet", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ws <- this$.ws;
  if (!is.null(ws))
    return(ws);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create weights files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create weights set. The data set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from data set");
  # Gets the WeightsSet Class object
  clazz <- getWeightsSetClass(this);
  ws <- clazz$fromDataSet(dataSet=ds, path=getPath(this), 
                                                     verbose=less(verbose));
  # make sure CDF is inherited
  setCdf(ws, getCdf(ds));
  verbose && exit(verbose);

  # Store in cache
  this$.ws <- ws;

  ws;
})

setMethodS3("getWeightsSetClass", "ProbeLevelModel", function(static, ...) {
  WeightsSet;
}, static=TRUE, private=TRUE)


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
setMethodS3("getFitFunction", "ProbeLevelModel", abstract=TRUE, static=TRUE, private=TRUE);



setMethodS3("getFitUnitFunction", "ProbeLevelModel", function(this, ...) {
  # Get the fit function for a single set of intensities
  fitfcn <- getFitFunction(this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the one for all blocks in a unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$probeModel == "pm-mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] - y[2,,];  # PM-MM
        fitfcn(y);
      })
    }
  } else if (this$probeModel == "min1(pm-mm)") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] - y[2,,];  # PM-MM
        y[y < 1] <- 1;       # min1(PM-MM)=min(PM-MM,1)
        fitfcn(y);
      })
    }
  } else if (this$probeModel == "pm+mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] + y[2,,];  # PM+MM
        fitfcn(y);
      })
    }
  } else {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        if (length(group) > 0) {
          y <- .subset2(group, 1); # Get intensities
        } else {
          y <- NULL;
        }
        fitfcn(y);
      })
    }
  }

  fitUnit;
}, private=TRUE)





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

  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ...);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);

  # Get the CEL intensities by units
  res <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && str(verbose, res[1]);
  verbose && exit(verbose);

  res;
}, private=TRUE)



setMethodS3("getCellIndices", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get what set of probes to read
  stratifyBy <- switch(this$probeModel, "pm"="pm", "mm"="mm", "pm-mm"="pmmm", "min1(pm-mm)"="pmmm", "pm+mm"="pmmm");

  # Get the CDF cell indices
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  verbose && enter(verbose, "Identifying CDF cell indices");
  verbose && cat(verbose, "Stratify by: ", stratifyBy);
  cells <- getCellIndices(cdf, ..., stratifyBy=stratifyBy, verbose=less(verbose));
  verbose && exit(verbose);
  
  cells;
}, private=TRUE)


###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
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
#  Returns an @integer @vector of unit indices.
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
setMethodS3("findUnitsTodo", "ProbeLevelModel", function(this, ...) {
  ces <- getChipEffectSet(this);
  findUnitsTodo(ces, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be fitted.
#     If @NULL, all units are considered.
#     If \code{remaining}, only non-fitted units are considered.
#   }
#   \item{...}{Arguments passed to @seemethod "readUnits".}
#   \item{force}{If @TRUE, already fitted units are re-fitted, and
#     cached data is re-read.}
#   \item{ram}{A @double indicating if more or less units should
#     be loaded into memory at the same time.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{moreUnits}{Deprected. Use \code{ram} instead.}
# }
#
# \value{
#  Returns an @integer @vector of indices of the units fitted, 
#  or @NULL if no units was (had to be) fitted.
# }
#
# \details{
#  All estimates are stored to file.
#
#  The non-array specific parameter estimates together with standard deviation
#  estimates and convergence information are stored in one file.
#
#  The parameter estimates specific to each array, typically "chip effects", 
#  are stored in array specific files.
#
#   Data set specific estimates [L = number of probes]:
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
setMethodS3("fit", "ProbeLevelModel", function(this, units="remaining", ..., force=FALSE, ram=moreUnits, verbose=FALSE, moreUnits=1) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cs <- getDataSet(this);
  cdf <- getCdf(cs);
  nbrOfArrays <- length(cs);

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

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(1e-3,Inf));

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

  # Get (and create if missing) the probe-affinity file (one per data set)
  paf <- getProbeAffinityFile(this, verbose=less(verbose));

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffectSet(this, verbose=less(verbose));

  # Number of units to load into memory and fit at the same time
  bytesPerChunk <- 100e6;       # 100Mb
  bytesPerUnitAndArray <- 500;  # Just a rough number; good enough?
  bytesPerUnit <- nbrOfArrays * bytesPerUnitAndArray;
  unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
  unitsPerChunk <- as.integer(max(unitsPerChunk,1));

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number units per chunk: ", unitsPerChunk);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Garbage collect
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  clearCache(cs);
  clearCache(paf);
  clearCache(ces);
  gc <- gc();
  verbose && print(verbose, gc);


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
    tRead <- processTime();
    y <- readUnits(this, units=units[uu], ..., force=force, cache=FALSE, 
                                                    verbose=less(verbose));
    timers$read <- timers$read + (processTime() - tRead);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting probe-level model");
    tFit <- processTime();
    fit <- lapply(y, FUN=fitUnit);
    timers$fit <- timers$fit + (processTime() - tFit);
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && str(verbose, fit[1]);
    verbose && exit(verbose);

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    verbose && print(verbose, gc);
    timers$gc <- timers$gc + (processTime() - tGc);

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

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    verbose && print(verbose, gc);
    timers$gc <- timers$gc + (processTime() - tGc);

    timers$total <- timers$total + (processTime() - tTotal);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # ETA
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose) {
      # Fraction left
      fLeft <- length(idxs) / nbrOfUnits;
      # Time this far
      lapTime <- processTime() - startTime;
      t <- Sys.time() - lapTime[3];
      printf(verbose, "Started: %s\n", format(t, "%Y%m%d %H:%M:%S"));
      # Estimated time left
      fDone <- 1-fLeft;
      timeLeft <- fLeft/fDone * lapTime;
      t <- timeLeft[3];
      printf(verbose, "Estimated time left: %.2fmin\n", t/60);
      # Estimate time to arrivale
      eta <- Sys.time() + t;
      printf(verbose, "ETA: %s\n", format(eta, "%Y%m%d %H:%M:%S"));
    }

    verbose && exit(verbose);
  } # while(length(idxs) > 0)


  totalTime <- processTime() - startTime;
  if (verbose) {
    nunits <- length(units);
    t <- totalTime[3];
    printf(verbose, "Total time for all units across all %d arrays: %.2fs == %.2fmin\n", nbrOfArrays, t, t/60);
    t <- totalTime[3]/nunits;
    printf(verbose, "Total time per unit across all %d arrays: %.2fs/unit\n", nbrOfArrays, t);
    t <- totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time per unit and array: %.3gms/unit & array\n", 1000*t);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", nbrOfUnits(cdf), t/60, t/3600);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits;
    printf(verbose, "Total time for complete data set: %.2fmin = %.2fh\n", t/60, t/3600);
    # Get distribution of what is spend where
    timers$write <- timers$writePaf + timers$writeCes;
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for writing chip-effects), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeCes"]/t["write"], t["gc"]);
  }

  invisible(units);
})




############################################################################
# HISTORY:
# 2007-04-02
# o Added support for the "pm+mm" probe model.
# 2007-02-29
# o BUG FIX: Probe-affinities was not save, resulting in all zeroes.
#   This was due to renaming getProbeAffinites() to getProbeAffinityFile().
# 2007-02-28
# o Added ETA to verbose output of fit() for the ProbeLevelModel.
# o Memory optimization: Further memory optimization by clearing the
#   cache of the 'cs', the 'paf', and the 'ces', before fitting.
# 2007-02-22
# o Added getChipEffectSet() and getProbeAffinityFile() to replace 
#   getChipEffects() and getProbeAffinites() in some future version.
# 2007-02-09
# o Added an additional garbage collection after fitting the PLM, but 
#   before storing parameter estimates.
# 2007-01-06
# o Now gc() memory information is outputted after each chunk.
# o Updated the formula for calculating the number of units per chunk in
#   fit(). Gives roughly the same number of units though.
# o Added probe model 'min1(PM-MM)' for modelling y = min(PM-MM,1), which
#   is how dChip v2006-12-14 is doing it.
# o Now ProbeLevelModel inherits directly from UnitModel.
# 2007-01-03
# o "Protected" several methods to simplify the interface for end users.
# o Added support from "PM-MM" probe models in addition to the default 
#   "PM only" model.
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
