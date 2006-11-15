###########################################################################/**
# @RdocClass AllelicCrosstalkCalibrator
#
# @title "The AllelicCrosstalkCalibrator class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the 
#  probe-level signals such that the signals from the two alleles are 
#  orthogonal.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#   \item{subversionTag}{A @character string appended to the version tag
#      of the input data set (with a period as a separator).}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
#   @include "../incl/AllelicCrosstalkCalibrator.Rex"
# }}
#
# @author
#*/###########################################################################
setConstructorS3("AllelicCrosstalkCalibrator", function(dataSet=NULL, subversionTag="", ...) {
  if (!is.null(dataSet)) {
    subversionTag <- Arguments$getCharacter(subversionTag);
    if (regexpr("[,.]", subversionTag) != -1)
      throw("A tag must not contain commas or periods: ", subversionTag);
  }

  extend(Object(), "AllelicCrosstalkCalibrator", 
    .subversionTag = subversionTag,
    inputDataSet = dataSet,
    "cached:outputDataSet" = NULL,
    .params = list(
    )
  )
})



setMethodS3("as.character", "AllelicCrosstalkCalibrator", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Input data set: %s", getName(ds)));
  tags <- paste(getTags(ds), collapse=", ");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Number of arrays: %d (%.2fMb)", 
                           nbrOfArrays(ds), getFileSize(ds)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Algorithm parameters: (%s)", params));
  s <- c(s, sprintf("Output path: %s", getPath(this)));
  s <- c(s, sprintf("Is done: %s", isDone(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("getSubversionTag", "AllelicCrosstalkCalibrator", function(this, ...) {
  this$.subversionTag;
})

setMethodS3("getVersionTag", "AllelicCrosstalkCalibrator", function(this, ...) {
  ds <- getInputDataSet(this);
  versionTag <- paste(getVersionTag(ds), getSubversionTag(this), sep=".");
  versionTag <- gsub("^[.]|[.]$", "", versionTag);
  versionTag;
})

setMethodS3("getName", "AllelicCrosstalkCalibrator", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getTags", "AllelicCrosstalkCalibrator", function(this, ...) {
  c(getVersionTag(this));
})

setMethodS3("getFullName", "AllelicCrosstalkCalibrator", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  fullname <- paste(name, tags, sep=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "AllelicCrosstalkCalibrator", function(this, ...) {
  "calibAllelicCT";
})


setMethodS3("getParametersAsString", "AllelicCrosstalkCalibrator", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, protected=TRUE)

setMethodS3("getParameters", "AllelicCrosstalkCalibrator", function(this, ...) {
  this$.params;
})

setMethodS3("getPath", "AllelicCrosstalkCalibrator", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  mkdirs(path);

  path;
})


###########################################################################/**
# @RdocMethod getInputDataSet
#
# @title "Gets the source data set"
#
# \description{
#  @get "title" that is to be normalized (or has been normalized).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSet", "AllelicCrosstalkCalibrator", function(this, ...) {
  this$inputDataSet;
})



###########################################################################/**
# @RdocMethod getOutputDataSet
#
# @title "Gets the normalized data set"
#
# \description{
#  @get "title", if normalized.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "AllelicCrosstalkCalibrator", function(this, ..., force=FALSE) {
  outputDataSet <- this$outputDataSet;
  if (force || is.null(outputDataSet)) {
    if (isDone(this)) {
      ds <- getInputDataSet(this);
      clazz <- Class$forName(class(ds)[1]);
      outputDataSet <- clazz$fromFiles(path=getPath(this));
      this$outputDataSet <- outputDataSet;
    }
  }
  outputDataSet;
})


setMethodS3("getOutputFiles", "AllelicCrosstalkCalibrator", function(this, ...) {
  outPath <- getPath(this);
  findFiles(pattern="[.](c|C)(e|E)(l|L)$", paths=outPath, firstOnly=FALSE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod isDone
#
# @title "Checks if the data set is calibrated or not"
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
#  Returns @TRUE if the data set is calibrated, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDone", "AllelicCrosstalkCalibrator", function(this, ...) {
  pathnames <- getOutputFiles(this);
  if (length(pathnames) == 0)
    return(FALSE);

  ds <- getInputDataSet(this);  
  if (length(pathnames) != nbrOfArrays(ds)) {
    throw("Number of output CEL files does not match the number of CEL files in the input dataset: ", length(pathnames), " != ", nbrOfArrays(ds));
  }
  
  return(TRUE);
})


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
# @examples "../incl/normalizeQuantile.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "AllelicCrosstalkCalibrator", function(this, ..., force=FALSE, verbose=FALSE) {
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

  # Already done?
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already calibrated");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input dataset
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- c(list(ds, path=outputPath, verbose=verbose), params);
  outputDataSet <- do.call("calibrateAllelicCrosstalk", args=args);

  # Update the output dataset
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})

############################################################################
# HISTORY:
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
