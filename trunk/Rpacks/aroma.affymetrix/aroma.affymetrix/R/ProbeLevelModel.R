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
#  "A [...] PLM is a model that is fit to probe-intensity data. 
#   More specifically, it is where we fit a model with probe level
#   and chip level parameters on a probeset by probeset basis",
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
#   \item{...}{Arguments passed to @see "AffymetrixUnitGroupsModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# 
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "ProbeAffinityFile" class.  
#   Use @seemethod "getProbeAffinities" to access these.
# }
#
# \seealso{
#   For more details on probe-level models, please see 
#   the \pkg{affyPLM} package.
# }
#*/###########################################################################
setConstructorS3("ProbeLevelModel", function(dataSet=NULL, path=filePath("modelPLM", getChipType(getCdf(dataSet))), model=c("pm"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'path':
  if (is.null(dataSet)) {
    # A work-around for the fact that getCdf(NULL) is not work.
    path=NULL;
  }

  extend(AffymetrixUnitGroupsModel(dataSet=dataSet, path=path, ...), "ProbeLevelModel",
    "cached:.paFile" = NULL,
    "cached:.chipFiles" = NULL,
    model = model
  )
}, abstract=TRUE)



###########################################################################/**
# @RdocMethod getProbeAffinityClass
#
# @title "Static method to get the ProbeAffinity Class object for this model"
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
  paFile <- this$.paFile;
  if (is.null(paFile)) {
    # Get the probe-affinity class object
    clazz <- getProbeAffinityClass(this);

    # Let the parameter object know about the CDF structure, because we 
    # might use a modified version of the one in the CEL header.
    cdf <- getCdf(this);
    pathname <- clazz$createFrom(cdf=cdf, path=getPath(this), dataSet=getDataSet(this));
    paFile <- newInstance(clazz, pathname, cdf=cdf, model=this$model);
    this$.paFile <- paFile;
  }
  paFile;
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
#  Returns a @character string.
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
setMethodS3("fit", "ProbeLevelModel", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the output directory etc
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The pathname of the CEL file where probe-affinity parameter estimates are stored
  paf <- getProbeAffinities(this);


  # For each of the data files, create a file to store the estimates in
  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  # Get (and create if missing) the chip files, which stores chip estimates
#  chipFiles <- getChipEffects(this);

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  stratifyBy <- switch(this$model, pm="pm");
  # IMPORTANT: This is the only place where the actually CDF structure
  # is being specified. The reading/writing is done using this structure!
  # This means that we can have getCellIndices() of AffymetrixCdfFile 
  # class to return a structure specific to our analysis.  For instance,
  # in CRLMM we estimate a model for each unit group, whereas in RLMM
  # and BRLMM, we estimate one model where the forward and reverse
  # strands have been joined. Thus, we only have to modify the 
  # getCellIndices() so it returns a structure that suits the.
  # model-fitting algorithm.  Pretty sweet actually :)  /HB 2006-08-24
  cdf0 <- getCellIndices(cdf, units=units, ..., stratifyBy=stratifyBy);
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading probe-affinity estimates");
  cel <- readUnits(paf, cdf=cdf0);
  verbose && str(verbose, cel, max.level=3, level=-5);
  verbose && exit(verbose);

  # Identify which of the requested units have *not* already been estimated
  if (force) {
    ntodo <- length(cel);
  } else {
    todo <- unlist(lapply(cel, FUN=function(unit) {
      any(isZero(.subset2(.subset2(unit, 1), "stdvs")));
    }), use.names=FALSE);
    ntodo <- sum(todo);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate model for set of units not already estimated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ntodo > 0) {
    fitfcn <- getFitFunction(this);

    if (force) {
      verbose && cat(verbose, "Forced estimation of ", ntodo, " units.");
      cdf <- cdf0;
    } else {
      verbose && cat(verbose, "Needs to estimate ", ntodo, " units.");
      cdf <- cdf0[todo];
    }

    # Get the CEL intensities by units
    verbose && enter(verbose, "Reading probe intensities");
    y <- getUnitIntensities(ds, units=cdf, ...);
    verbose && exit(verbose);
  
    # Apply the Li-Wong model
    verbose && enter(verbose, "Fitting probeset model");
    fit <- lapply(y, FUN=function(unit) {
      groups <- lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1);
        f <- fitfcn(y);
      })
      groups;
    })
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing parameter estimates");
    fit <- updateUnits(paf, cdf=cdf, data=fit);
    verbose && exit(verbose);
    fit <- decode(paf, fit);

    verbose && enter(verbose, "Mixing parameter estimates");
    if (force) {
      cel <- fit;
    } else {
      cel[todo] <- fit;
    }
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "No need to fit anything. Found estimates for all units requested.");
  } # if (length(cdf) > 0)

  cel;
})

setMethodS3("getChipEffects", "ProbeLevelModel", function(this, ..., verbose=FALSE) {
  chipFiles <- this$.chipFiles;
  if (!is.null(chipFiles))
    return(chipFiles);

  # For each of the data files, create a file to store the estimates in
  path <- getPath(this);
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  unitSizes <- getUnitSizes(cdf);
  unitOffsets <- cumsum(unitSizes) - unitSizes[1] + 1;
  n <- sum(unitSizes);

  chipFiles <- list();
  for (kk in seq(ds)) {
    df <- ds[[kk]];
    filename <- sprintf("%s-liwong.apd", getName(df));
    filename <- filePath(path, filename);
    if (!isFile(filename)) {
      X <- FileFloatVector(filename, length=n);
      X <- FileFloatVector(length=n, appendTo=X);
      X <- FileByteVector(length=n, appendTo=X);
      close(X);
      rm(X);
    }
    set <- AbstractFileArray$fromFile(filename);
    # We have to close the parameter files, because Windows can only
    # handle ~512 open connections, and we might have more arrays.
    # Instead we have to open and close the connections each time we
    # read data.
    close(set[[1]]);
    chipFiles[[kk]] <- set;
  } # for (kk in ...)

  this$.chipFiles <- chipFiles;

  chipFiles;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-08-25
# o Created from AffymetrixLiWongModel.  So much is in common with the RMA
#   model so a lot can be reused if done cleverly.
############################################################################
