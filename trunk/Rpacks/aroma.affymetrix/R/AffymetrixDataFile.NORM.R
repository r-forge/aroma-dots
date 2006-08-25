###########################################################################/**
# @set "class=AffymetrixDataFile"
# @RdocMethod normalizeQuantile
#
# @title "Normalizes the sample to a target empirical distribution"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{outPath}{The path where to save the normalized data files.}
#   \item{xTarget}{A @numeric @vector.  The empirical distribution
#     to which all arrays should be normalized to.}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
#   \item{overwrite}{If @TRUE, already normalized arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{format}{The output file format.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the normalized @see "AffymetrixDataFile" object.
# }
#
# @author
#
# \seealso{
#   @see "aroma.light::normalizeQuantile.numeric"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("normalizeQuantile", "AffymetrixDataFile", function(this, outPath=file.path("norm", getChipType(this)), xTarget, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, format=getOption("aroma.affymetrix/format"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  outPath <- Arguments$getWritablePathname(outPath);
  if (identical(getPath(this), outPath)) {
    throw("Cannot not normalize data file. Argument 'outPath' refers to the same path as the path of the data file to be normalized: ", outPath);
  }
  mkdirs(outPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'format':
  if (is.null(format)) {
    format <- "apd";
  } else if (format == "apd") {
  } else if (format == "cel") {
    if (getFileType(this) != "cel") {
      throw("Cannot not normalize data.  Argument 'format' is 'cel', but currently only normalized data originating from CEL files can be written as CEL files: ", getFileType(this));
    }
  } else {
    throw("Unknown value of argument 'format': ", format);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- paste(getName(this), format, sep=".");
  pathname <- Arguments$getWritablePathname(filename, path=outPath, 
                                         mustNotExist=(!overwrite && !skip));

  # Already normalized?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Normalized data file already exists: ", pathname);
    return(AffymetrixDataFile$fromFile(pathname));
  }

  # Get probe signals
  x <- getProbeIntensities(this, ..., verbose=verbose);

  # Identify the subset of probes to be updated
  subsetToUpdate <- getProbes(this, probes=subsetToUpdate, type=typesToUpdate, verbose=verbose);

  # Normalize intensities
  verbose && enter(verbose, "Normalizing to empirical target distribution");
  x[subsetToUpdate] <- normalizeQuantile(x[subsetToUpdate], xTarget=xTarget);
  rm(subsetToUpdate);
  verbose && exit(verbose);

  # Write normalized data to file
  verbose && enter(verbose, "Writing normalized probe signals");
  if (format == "apd") {
    chipType <- getChipType(this);
    apdMap <- getApdMap(this);
    if (is.null(apdMap)) {
      mapType <- NULL;
      writeMap <- NULL;
    } else {
      mapType <- getMapType(apdMap);
      writeMap <- getWriteMap(apdMap);
    }
    rm(apdMap);
    if (isFile(pathname) && overwrite) {
      file.remove(pathname);
    }
    writeApd(pathname, data=x, chipType=chipType, mapType=mapType, writeMap=writeMap);
  } else {
    # Copy CEL file and update the copy
    copyCel(from=getPathname(this), to=pathname, overwrite=overwrite);
    updateCel(pathname, intensities=x);
  }
  rm(x);
  verbose && exit(verbose);

  # Return normalized data file object
  AffymetrixDataFile$fromFile(pathname);
})



###########################################################################/**
# @RdocMethod fitQuantileNormFcn
#
# @title "Fits quantile normalization functions for the arrays in the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{yTarget}{The target probe signals.}
#   \item{subset}{An optional @numeric @vector specifying a subset of probe
#      indices used to fit the normalization function.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a normalization @function.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitQuantileNormFcn", "AffymetrixDataFile", function(this, yTarget, subset=NULL, ..., controlParams=list(spar=NULL, nknots=1024), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'yTarget':
  if (!identical(attr(yTarget, "isSorted"), TRUE)) {
    # Sort target signals
    verbose && enter(verbose, "Sorting ", length(yTarget), " target signals");
    yTarget <- sort(yTarget);
    attr(yTarget, "isSorted") <- TRUE;
    verbose && exit(verbose);
  }

  # Argument 'controlParams':
  spar <- controlParams$spar;
  nknots <- controlParams$nknots;

  verbose && enter(verbose, "Fitting (quantile) normalization function");

  # Read the probe intensities
  y <- readIntensities(this, probes=subset, verbose=verbose);

  # Sort signals
  verbose && enter(verbose, "Sorting probe signals");
  y <- sort(y);
  verbose && exit(verbose);

  # Fit normalization function
  verbose && enter(verbose, "Fitting smooth spline");
  ok <- !is.na(yTarget) & !is.na(y);
  sp <- smooth.spline(x=y[ok], y=yTarget[ok], spar=spar, nknots=nknots, 
                                                          keep.data=FALSE);
  verbose && exit(verbose);

  # Create a minimal 'smooth.spline' object for prediction.
  # Note: You can not do predict(sp$fit, ...) because then there
  # will be a problem with Recall().  See r-devel on 2006-04-05.
  fit <- structure(list(fit=sp$fit), class=class(sp))

  env <- new.env(parent=baseenv());
  assign("fit", fit, envir=env);
  # Create transformation function
  fcn <- function(x, ...) {
    stats::predict(fit, x, ...)$y;
  }
  environment(fcn) <- env;

  verbose && exit(verbose);

  fcn;
}, private=TRUE) #fitQuantileNormFcn()




setMethodS3("transformAffine", "AffymetrixDataFile", function(this, outPath=file.path("transformed", getChipType(this)), offset=0, scale=1, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, format=getOption("aroma.affymetrix/format"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  outPath <- Arguments$getWritablePathname(outPath);
  if (identical(getPath(this), outPath)) {
    throw("Cannot not transform data. Argument 'outPath' refers to the same path as the path of the data file to be transformed: ", outPath);
  }
  mkdirs(outPath);

  # Argument 'offset':
  offset <- Arguments$getDouble(offset);

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'format':
  if (is.null(format)) {
    format <- "apd";
  } else if (format == "apd") {
  } else if (format == "cel") {
    if (getFileType(this) != "cel") {
      throw("Cannot not transform data.  Argument 'format' is 'cel', but currently only transform data originating from CEL files can be written as CEL files: ", getFileType(this));
    }
  } else {
    throw("Unknown value of argument 'format': ", format);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- paste(getName(this), format, sep=".");
  pathname <- Arguments$getWritablePathname(filename, path=outPath, 
                                         mustNotExist=(!overwrite && !skip));

  # Already shifted?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Transformed data file already exists: ", pathname);
    return(AffymetrixDataFile$fromFile(pathname));
  }

  # Get probe signals
  x <- getProbeIntensities(this, ..., verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subsetToUpdate <- getProbes(this, probes=subsetToUpdate, type=typesToUpdate, verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift intensities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Transforming probe intensities by (offset,scale)=(%.1f,%.2f) ", offset, scale));
  x[subsetToUpdate] <- offset + scale*x[subsetToUpdate];
  rm(subsetToUpdate);
  verbose && exit(verbose);

  # Write normalized data to file
  verbose && enter(verbose, "Writing transformed probe signals");
  if (format == "apd") {
    chipType <- getChipType(this);
    apdMap <- getApdMap(this);
    if (is.null(apdMap)) {
      mapType <- NULL;
      writeMap <- NULL;
    } else {
      mapType <- getMapType(apdMap);
      writeMap <- getWriteMap(apdMap);
    }
    rm(apdMap);
    if (isFile(pathname) && overwrite) {
      file.remove(pathname);
    }
    writeApd(pathname, data=x, chipType=chipType, mapType=mapType, writeMap=writeMap);
  } else {
    # Copy CEL file and update the copy
    copyCel(from=getPathname(this), to=pathname, overwrite=overwrite);
    updateCel(pathname, intensities=x);
  }
  rm(x);
  verbose && exit(verbose);

  # Return transformed data file object
  AffymetrixDataFile$fromFile(pathname);
}) # transformAffine()

############################################################################
# HISTORY:
# 2006-07-28
# o Added transformAffine().
# 2006-07-21
# o Added more verbose output for normalizeQuantile().
# o BUG FIX: typo for normalizeQuantile(..., format="cel").
# 2006-07-08
# o Added argument 'format' to normalizeQuantile().
# o Added support to write normalized data to CEL files.  Currently we do
#   this by copying the existing CEL file and updating that.
# 2006-05-15
# o Created from AffymetrixDataFile.R.
# 2006-04-05
# o BUG FIX:  fitQuantileNormFcn() returned a function that when called with
#   x values forcing extrapolation, error "Error in Recall(object, xrange) :
#   couldn't find function "predict.smooth.spline.fit" would be thrown.
#   This is because you cannot do predict(sp$fit, ...) but only 
#   predict(sp, ...).  Why I don't really know; probably about namespaces.
# 2006-03-18
# o Added argument 'subset' to fitQuantileNormFcn().
# 2006-03-03
# o When creating transformation function in, say, fitQuantileNormFcn(), it
#   is important to create an empty environment for the function otherwise
#   all arguments in the calling function is included too.
# o Added writeApd().  For now it can only write 'intensities'.
############################################################################
