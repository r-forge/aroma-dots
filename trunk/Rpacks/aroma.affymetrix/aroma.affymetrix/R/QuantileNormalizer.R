###########################################################################/**
# @RdocClass QuantileNormalizer
#
# @title "The QuantileNormalizer class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization function that transforms the 
#  probe-level signals towards the same empirical distribution.
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
#   @include "../incl/QuantileNormalizer.Rex"
# }}
#
# @author
#*/###########################################################################
setConstructorS3("QuantileNormalizer", function(dataSet=NULL, subversionTag="", subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, ...) {
  if (!is.null(dataSet)) {
    subversionTag <- Arguments$getCharacter(subversionTag);
    if (regexpr("[,.]", subversionTag) != -1)
      throw("A tag must not contain commas or periods: ", subversionTag);
  }

  extend(Object(), "QuantileNormalizer", 
    .subversionTag = subversionTag,
    inputDataSet = dataSet,
    "cached:outputDataSet" = NULL,
    .params = list(
       subsetToUpdate = subsetToUpdate,
       typesToUpdate = typesToUpdate,
       subsetToAvg = subsetToAvg,
       typesToAvg = typesToAvg,
       .targetDistribution = targetDistribution
    )
  )
})



setMethodS3("as.character", "QuantileNormalizer", function(this, ...) {
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


setMethodS3("getSubversionTag", "QuantileNormalizer", function(this, ...) {
  this$.subversionTag;
})

setMethodS3("getVersionTag", "QuantileNormalizer", function(this, ...) {
  ds <- getInputDataSet(this);
  versionTag <- paste(getVersionTag(ds), getSubversionTag(this), sep=".");
  versionTag <- gsub("^[.]|[.]$", "", versionTag);
  versionTag;
})

setMethodS3("getName", "QuantileNormalizer", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getTags", "QuantileNormalizer", function(this, ...) {
  c(getVersionTag(this));
})

setMethodS3("getFullName", "QuantileNormalizer", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  fullname <- paste(name, tags, sep=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "QuantileNormalizer", function(this, ...) {
  "normQuantile";
})


setMethodS3("getParametersAsString", "QuantileNormalizer", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, protected=TRUE)

setMethodS3("getParameters", "QuantileNormalizer", function(this, ...) {
  this$.params;
})

setMethodS3("getOutputIdentifier", "QuantileNormalizer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calculating the output identifier");

  verbose && enter(verbose, "Retrieving the identifier for input data set");
  ds <- getInputDataSet(this);
  inputId <- getIdentifier(ds);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating the identifier for parameters");
  paramId <- this$.paramId;
  params <- getParameters(this);
  params$.targetDistribution <- attr(params$.targetDistribution, "identifier");
  paramId <- digest(list(params));
  this$.paramId <- paramId;
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating the joint identifier");
  id <- digest(list(inputId, paramId));
  verbose && exit(verbose);
 
  verbose && exit(verbose);

  id;
}, protected=TRUE)



setMethodS3("getPath", "QuantileNormalizer", function(this, ...) {
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
setMethodS3("getInputDataSet", "QuantileNormalizer", function(this, ...) {
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
setMethodS3("getOutputDataSet", "QuantileNormalizer", function(this, ..., force=FALSE) {
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


setMethodS3("getOutputFiles", "QuantileNormalizer", function(this, ...) {
  outPath <- getPath(this);
  findFiles(pattern="[.](c|C)(e|E)(l|L)$", paths=outPath, firstOnly=FALSE);
}, protected=TRUE)


setMethodS3("getTargetDistribution", "QuantileNormalizer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting target distribution");

  yTarget <- this$.params$.targetDistribution;
  if (force || is.null(yTarget)) {
    pathname <- getTargetDistributionPathname(this, verbose=less(verbose));
    verbose && print(verbose, pathname);

    if (isFile(pathname)) {
      verbose && enter(verbose, "Reading saved distribution: ", pathname);
      yTarget <- readApd(pathname)$quantiles;
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Calculating");
      yTarget <- calculateTargetDistribution(this, verbose=less(verbose));
      verbose && exit(verbose);
    }
    attr(yTarget, "identifier") <- getTargetDistributionIdentifier(this);
    this$.params$.targetDistribution <- yTarget;
  } else {
    verbose && cat(verbose, "Was cached in-memory.");
  }

  verbose && exit(verbose);

  yTarget;
})


setMethodS3("getTargetDistributionIdentifier", "QuantileNormalizer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting identifier for target distribution");

  ds <- getInputDataSet(this);
  params <- getParameters(this);
  # Get the parameters used for averaging
  id <- digest(list(
    identifier=getIdentifier(ds), 
    indices=params$subsetToAvg, 
    types=params$typesToAvg
  ));

  verbose && exit(verbose);

  id;
}, protected=TRUE)

setMethodS3("getTargetDistributionPathname", "QuantileNormalizer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting pathname for target distribution");

  id <- getTargetDistributionIdentifier(this, verbose=less(verbose));
  ds <- getInputDataSet(this);
  filename <- sprintf(".averageQuantile-%s.apq", id);
  path <- getPath(ds);
  pathname <- filePath(path, filename, expandLinks="any");

  verbose && exit(verbose);

  pathname;
}, protected=TRUE)

setMethodS3("calculateTargetDistribution", "QuantileNormalizer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating target distribution");
  verbose && cat(verbose, "Method: average empirical distribution");

  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  params <- getParameters(this);
  probes <- identifyCells(cdf, indices=params$subsetToAvg, 
                         types=params$typesToAvg, verbose=less(verbose));
  verbose && cat(verbose, "Using ", length(probes), " probes");
  verbose && cat(verbose, "Calculating target distribution from the ", length(ds), " arrays in the input data set");

  # Calculate the average quantile
  yTarget <- averageQuantile(ds, probes=probes, verbose=less(verbose));

  # Write the result to file
  pathname <- getTargetDistributionPathname(this, verbose=less(verbose));
  verbose && cat(verbose, "Saving distribution: ", pathname);
  writeApd(pathname, data=yTarget, name="quantiles");

  verbose && exit(verbose);

  invisible(yTarget);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod isDone
#
# @title "Checks if the data set is normalized or not"
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
#  Returns @TRUE if the data set is normalized, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDone", "QuantileNormalizer", function(this, ...) {
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
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
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
setMethodS3("process", "QuantileNormalizer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Quantile normalizing data set");

  # Already done?
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve/calculate the target distribution
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getTargetDistribution(this, verbose=less(verbose));

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
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  names(params) <- gsub(".targetDistribution", "xTarget", names(params));
  args <- c(list(ds, path=outputPath, verbose=verbose), params);
  outputDataSet <- do.call("normalizeQuantile", args=args);

  # Update the output dataset
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})

############################################################################
# HISTORY:
# 2006-10-30
# o Created.
############################################################################
