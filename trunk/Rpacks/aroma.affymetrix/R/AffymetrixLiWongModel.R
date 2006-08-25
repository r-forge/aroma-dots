###########################################################################/**
# @RdocClass AffymetrixLiWongModel
#
# @title "The AffymetrixLiWongModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Li \& Wong (2001) model.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixUnitGroupsModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "LiWongProbeAffinityFile" class.  
#   Use @seemethod "getProbeAffinities" to access these.
# }
#
# \references{
#   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
#   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
# }
#*/###########################################################################
setConstructorS3("AffymetrixLiWongModel", function(dataSet=NULL, path=filePath("modelLiWong", getChipType(getCdf(dataSet))), model=c("pm"), ...) {
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

  extend(AffymetrixUnitGroupsModel(dataSet=dataSet, path=path, ...), "AffymetrixLiWongModel",
    "cached:.paFile" = NULL,
    "cached:.chipFiles" = NULL,
    model = model
  )
})



setMethodS3("getProbeAffinities", "AffymetrixLiWongModel", function(this, ...) {
  paFile <- this$.paFile;
  if (is.null(paFile)) {
    # Let the parameter object know about the CDF structure, because we 
    # might use a modified version of the one in the CEL header.
    cdf <- getCdf(this);
    pathname <- LiWongProbeAffinityFile$createFrom(cdf=cdf, path=getPath(this), dataSet=getDataSet(this));
    paFile <- LiWongProbeAffinityFile(pathname, cdf=cdf, model=this$model);
    this$.paFile <- paFile;
  }
  paFile;
})


setMethodS3("getChipEffects", "AffymetrixLiWongModel", function(this, ..., verbose=FALSE) {
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
# @keyword programming
#*/###########################################################################
setMethodS3("fit", "AffymetrixLiWongModel", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  require(affy) || throw("Package 'affy' not loaded.");

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
      any(isZero(.subset2(.subset2(unit, 1), "nbrOfIterations")));
    }), use.names=FALSE);
    ntodo <- sum(todo);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate model for set of units not already estimated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ntodo > 0) {
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
        y <- t(.subset2(group, 1));
        f <- fit.li.wong(y);
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



setMethodS3("setup", "AffymetrixLiWongModel", function(this, path=file.path("liwong", getChipType(this)), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePathname(path, mustExist=FALSE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the output directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the file to store probe-affinities etc
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- Arguments$getWritablePathname("probe-affinities.CEL", path=path);
  if (!isFile(pathname)) {
    verbose && enter(verbose, "Creating file for probe-affinity estimates.");
    file.copy(getPathname(ds[[1]]), pathname);
    # Clear the file
    updateCel(pathname, pixels=rep(0:0, length=nbrOfCells(getCdf(this))));
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the files to store array-specific estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  for (kk in seq(ds)) {
    df <- ds[[kk]];
    filename <- sprintf("%s.apd", getName(df));
    pathname <- file.path(path, filename);
    if (!isFile(pathname)) {
      verbose && enter(verbose, "Creating file for array-specific estimates.");
      verbose && exit(verbose);
    }
  }
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added getProbeAffinities() and the corrsponding cached fields.
# o Now fit() does not re-read data just updated.
# 2006-08-19
# o After all the bug fixes in updateCel() I think this function finally
#   works correctly.
# 2006-08-17
# o Created.
############################################################################
