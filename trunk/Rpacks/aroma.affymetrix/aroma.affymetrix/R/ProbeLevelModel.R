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

  extend(AffymetrixUnitGroupsModel(..., name=name), "ProbeLevelModel",
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

setMethodS3("getChipEffectClass", "ProbeLevelModel", function(static, ...) {
  ChipEffectFile;
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
setMethodS3("fit", "ProbeLevelModel", function(this, units="remaining", ..., transform=NULL, postCdfTransform=NULL, unitsPerChunk=moreUnits*100000/length(getDataSet(this)), moreUnits=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  # Get (and create if missing) the probe-affinity file (one per dataset)
  paf <- getProbeAffinities(this);

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffects(this);


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
    units <- findUnitsTodo(paf);
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'transform':
  if (is.null(transform)) {
  } else if (is.function(transform)) {
  } else {
    throw("Argument 'transform' is not a function: ", mode(transform));
  }

  # Argument 'postCdfTransform':
  if (is.null(postCdfTransform)) {
  } else if (is.function(postCdfTransform)) {
  } else {
    throw("Argument 'postCdfTransform' is not a function: ", mode(postCdfTransform));
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
      units <- findUnitsTodo(paf, units=units, transform=postCdfTransform);
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
  # Get what set of probes to read
  stratifyBy <- switch(this$model, pm="pm");

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

    # Get the CDF cell indices
    verbose && enter(verbose, "Identifying CDF cell indices");
    cdfUnits <- getCellIndices(cdf, units=units[uu], ..., stratifyBy=stratifyBy);
    verbose && exit(verbose);

    # Get the CEL intensities by units
    verbose && enter(verbose, "Reading probe intensities");
    y <- getUnitIntensities(ds, units=cdfUnits, ...);
    verbose && exit(verbose);
  
    # Fit the model
    verbose && enter(verbose, "Fitting unit-group model");
    fit <- lapply(y, FUN=function(unit) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1);
        if (!is.null(transform))
          y <- transform(y);
        fitfcn(y);
      })
    })
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing chip-effect estimates");
    firstCells <- getFirstCellIndices(paf, units=units[uu]);
    if (!is.null(postCdfTransform)) {
      firstCells <- postCdfTransform(firstCells);
    }
    updateUnits(ces, cdf=firstCells, data=fit, verbose=verbose);
    firstCells <- NULL; # Not needed anymore
    verbose && exit(verbose);

    # Store the probe-affinities *last* because we use these to 
    # identify which units are already modelled or not.
    verbose && enter(verbose, "Storing probe-affinity estimates");
    if (!is.null(postCdfTransform)) {
      cdfUnits <- postCdfTransform(cdfUnits);
    }
    updateUnits(paf, cdf=cdfUnits, data=fit);
    cdfUnits <- fit <- NULL; # Not needed anymore
    verbose && exit(verbose);

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
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  chipFiles <- this$.chipFiles;
  if (!is.null(chipFiles))
    return(chipFiles);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe-affinity class object
  clazz <- getChipEffectClass(this);

  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  cdf <- getCdf(this);
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect file. The CEL set is empty.");
  
  # For each CEL file, create a chip-effect file
  df <- as.list(ds)[[1]];

  template <- NULL;  
  chipFiles <- list();
  nFiles <- length(ds);
  for (kk in seq(ds)) {
    df <- as.list(ds)[[kk]];
    name <- getName(df);
    verbose && enter(verbose, sprintf("Creating chip-effect file #%d of %d: %s", kk, nFiles, name));
    filename <- sprintf("%s-chipEffect.CEL", name);
    pathname <- Arguments$getWritablePathname(filename, path=getPath(this));

    # Use template to create missing files
    if (!isFile(pathname)) {
      # Create an empty template CEL file for quick copy of the others
      if (is.null(template)) {
        verbose && enter(verbose, "Creating template");
        template <- createFrom(df, filename=".template.CEL", path=getPath(this), 
                                                               verbose=verbose);
        template <- getPathname(template);
        on.exit(file.remove(template));
        verbose && exit(verbose);
      }
  
      verbose && enter(verbose, "Copying template");
      file.copy(template, pathname);
      verbose && exit(verbose);
    } else {
      verbose && cat(verbose, "File already exists.");
    }

    # Create a chip-effect object
    chipFile <- newInstance(clazz, pathname, cdf=cdf, model=this$model);

    chipFiles <- c(chipFiles, list(chipFile));

    verbose && exit(verbose);
  }

  chipFiles <- ChipEffectSet(chipFiles);
  this$.chipFiles <- chipFiles;

  chipFiles;
})


############################################################################
# HISTORY:
# 2006-08-26
# o Added new getChipEffects().
# 2006-08-25
# o Created from AffymetrixLiWongModel.  So much is in common with the RMA
#   model so a lot can be reused if done cleverly.
############################################################################
