###########################################################################/**
# @RdocClass QuantileNormalization
#
# @title "The QuantileNormalization class"
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
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{targetDistribution}{A @numeric @vector.  The empirical 
#     distribution to which all arrays should be normalized to.}
#   \item{subsetToAvg}{The probes to calculate average empirical
#     distribution over.  If a single @numeric in (0,1), then this
#     fraction of all probes will be used.  
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the 
#     average empirical distribution.  
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch 
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
#   @include "../incl/QuantileNormalization.Rex"
# }}
#
# @author
#*/###########################################################################
setConstructorS3("QuantileNormalization", function(..., subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate) {
  extend(ProbeLevelTransform(...), "QuantileNormalization", 
    .subsetToUpdate = subsetToUpdate,
    .typesToUpdate = typesToUpdate,
    .targetDistribution = targetDistribution,
    .subsetToAvg = subsetToAvg,
    .typesToAvg = typesToAvg
  )
})


setMethodS3("getSubsetToUpdate", "QuantileNormalization", function(this, ...) {
  this$.subsetToUpdate;
}, private=TRUE)


setMethodS3("getParameters", "QuantileNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    subsetToUpdate = this$.subsetToUpdate,
    typesToUpdate = this$.typesToUpdate,
    subsetToAvg = this$.subsetToAvg,
    typesToAvg = this$.typesToAvg,
    .targetDistribution = this$.targetDistribution
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)



setMethodS3("getTargetDistribution", "QuantileNormalization", function(this, sort=TRUE, ..., force=FALSE, verbose=FALSE) {
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

  yTarget <- this$.targetDistribution;
  if (inherits(yTarget, "AffymetrixCelFile")) {
    df <- yTarget;
    verbose && enter(verbose, "Reading distribution from baseline array");
    verbose && cat(verbose, "Array: ", getFullName(df));
    yTarget <- getData(df, field="intensities")$intensities;
    this$.targetDistribution <- yTarget;
    verbose && exit(verbose);
  } else if (force || is.null(yTarget)) {
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
    this$.targetDistribution <- yTarget;
  } else {
    verbose && cat(verbose, "Was specified or cached in-memory.");
  }

  if (sort)
    yTarget <- sort(yTarget);
  verbose && str(verbose, yTarget);

  verbose && exit(verbose);

  yTarget;
}, private=TRUE)


setMethodS3("getTargetDistributionIdentifier", "QuantileNormalization", function(this, ..., verbose=FALSE) {
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
  cdf <- getCdf(ds);
  params <- getParameters(this);
  # Get the parameters used for averaging
  indices <- params$subsetToAvg;

  # Speed up for digest() in case there are many indices
  nbrOfCells <- nbrOfCells(cdf);
  if (length(indices) > nbrOfCells/2) {
    indices <- -setdiff(1:nbrOfCells, indices);
  }

  key <- list(
    identifier=getIdentifier(ds), 
    indices=indices,
    types=params$typesToAvg
  );
  id <- digest(key);
  verbose && exit(verbose);

  id;
}, private=TRUE)

setMethodS3("getTargetDistributionPathname", "QuantileNormalization", function(this, ..., verbose=FALSE) {
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

  ds <- getInputDataSet(this);
  path <- getPath(ds);
  id <- getTargetDistributionIdentifier(this, verbose=less(verbose));
  filename <- sprintf(".averageQuantile-%s.apq", id);
  pathname <- filePath(path, filename, expandLinks="any");

  verbose && exit(verbose);

  pathname;
}, private=TRUE)


setMethodS3("calculateTargetDistribution", "QuantileNormalization", function(this, ..., verbose=FALSE) {
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

  pathname <- getTargetDistributionPathname(this, verbose=less(verbose));

  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  params <- getParameters(this);
  probes <- identifyCells(cdf, indices=params$subsetToAvg, 
                         types=params$typesToAvg, verbose=less(verbose));
  rm(params);
  verbose && cat(verbose, "Using ", length(probes), " probes");
  verbose && cat(verbose, "Calculating target distribution from the ", length(ds), " arrays in the input data set");

  # Calculate the average quantile
  yTarget <- averageQuantile(ds, probes=probes, verbose=less(verbose));
  rm(probes);

  # Write the result to file
  verbose && cat(verbose, "Saving distribution: ", pathname);
  writeApd(pathname, data=yTarget, name="quantiles");

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  invisible(yTarget);
}, private=TRUE)




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
setMethodS3("process", "QuantileNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve/calculate the target distribution
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving target distribution");
  getTargetDistribution(this, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters (including the target distribution)
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing data towards target distribution");
  names(params) <- gsub(".targetDistribution", "xTarget", names(params));
  args <- c(list(ds, path=outputPath, verbose=verbose), params);

  # Garbage collect
  rm(params); gc();

  outputDataSet <- do.call("normalizeQuantile", args=args);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})

############################################################################
# HISTORY:
# 2007-02-04
# o Now QuantileNormalization() takes an AffymetrixCelFile as a target
#   distribution too, cf argument 'targetDistribution'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessor class.
# o Now this pre-processor output results to probeData/.
# o Renamed from QuantileNormalizer.
# 2006-11-18
# o Removed version and subversion tags, and related functions. 
#   Now getTags() returns the tags of the input data set plus any tags 
#   of this instance.
# 2006-10-30
# o Created.
############################################################################
