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
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
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
setConstructorS3("AllelicCrosstalkCalibrator", function(dataSet=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", class(dataSet));
  }

  # Argument 'tags':
  if (!is.null(dataSet)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "ACT";
  }


  extend(Object(), "AllelicCrosstalkCalibrator", 
    .tags = tags,
    inputDataSet = dataSet,
    "cached:outputDataSet" = NULL,
    .params = list(
    )
  )
})



setMethodS3("as.character", "AllelicCrosstalkCalibrator", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  tags <- paste(getTags(ds), collapse=", ");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
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


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "AllelicCrosstalkCalibrator", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this normalizer.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "AllelicCrosstalkCalibrator", function(this, ...) {
  ds <- getInputDataSet(this);
  c(getTags(ds), this$.tags);
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "AllelicCrosstalkCalibrator", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "AllelicCrosstalkCalibrator", function(this, ...) {
  "calibAllelicCT";
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
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
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

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
  # Get input dataset
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
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
# 2006-11-18
# o Removed version and subversion tags, and related functions. 
#   Now getTags() returns the tags of the input data set plus any tags 
#   of this instance.
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
