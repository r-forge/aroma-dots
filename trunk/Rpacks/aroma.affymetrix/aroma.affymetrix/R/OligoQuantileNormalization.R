###########################################################################/**
# @RdocClass OligoQuantileNormalization
#
# @title "The OligoQuantileNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a special @see "QuantileNormalization" that produces
#  identical output as the quantile normalization in the \pkg{oligo} package.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#       @see "QuantileNormalization".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("OligoQuantileNormalization", function(...) {
  extend(QuantileNormalization(..., typesToUpdate="pm"), "OligoQuantileNormalization");
})


setMethodS3("getSubsetToUpdate", "OligoQuantileNormalization", function(this, ..., verbose=FALSE) {
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);
  subset <- which(isPm(pd));
  cells <- getCellIndices(pd, subset=subset);
  cells;
})


setMethodS3("getTargetDistribution", "OligoQuantileNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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
  if (force || is.null(yTarget)) {
    ds <- getInputDataSet(this);
    cdf <- getCdf(ds);
    pd <- PlatformDesign(cdf);
    yTarget <- getReferenceQuantiles(pd);
    this$.targetDistribution <- yTarget;
  } else {
    verbose && cat(verbose, "Was specified or cached in-memory.");
  }

  verbose && exit(verbose);

  yTarget;
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "OligoQuantileNormalization", function(this, ..., force=FALSE, skip=TRUE, verbose=FALSE) {
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
#  if (!force && isDone(this)) {
#    verbose && cat(verbose, "Already normalized");
#    verbose && exit(verbose);
#    outputDataSet <- getOutputDataSet(this);
#    return(invisible(outputDataSet));
#  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Retrieve/calculate the target distribution
  yTarget <- getTargetDistribution(this, verbose=less(verbose));

  # Get algorithm parameters
  subsetToUpdate <- getSubsetToUpdate(this);

  if (length(yTarget) != length(subsetToUpdate)) {
    throw(sprintf("Error in platform-design package '%s'. The number of PM values in reference distribution does not match the number of PM features: %s != %s", getPackageName(pd), length(yTarget), length(subsetToUpdate)));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing ", nbrOfArrays(ds), " arrays");
  dataFiles <- list();
  for (kk in seq(ds)) {
    verbose && enter(verbose, "Array #", kk);
    df <- getFile(ds, kk);
    verbose && print(verbose, df);

    filename <- basename(getPathname(df));
    pathname <- Arguments$getWritablePathname(filename, path=outputPath);
  
    # Already normalized?
    if (isFile(pathname) && skip) {
      verbose && cat(verbose, "Normalized data file already exists: ", pathname);
      # CDF inheritance
      dataFiles[[kk]] <- fromFile(df, pathname);
      verbose && exit(verbose);
      next;
    }
  
    # Get all probe signals
    verbose && enter(verbose, "Reading probe intensities");
    x <- getData(df, indices=subsetToUpdate, fields="intensities", verbose=less(verbose,2));
    x <- x$intensities;
    verbose && str(verbose, x);
    verbose && exit(verbose);
  
    x <- oligo::normalizeToSample(as.matrix(x), yTarget)[,1];

    # Write normalized data to file
    verbose && enter(verbose, "Writing normalized probe signals");
    # Copy CEL file and update the copy
    verbose && enter(verbose, "Copying source CEL file");
    copyCel(from=getPathname(df), to=pathname, overwrite=!skip);
    verbose && exit(verbose);
    verbose && enter(verbose, "Writing normalized intensities");
    updateCel(pathname, indices=subsetToUpdate, intensities=x);
    rm(x);
    verbose && exit(verbose);
    verbose && exit(verbose);
  
    # Return new normalized data file object
    dataFiles[[kk]] <- fromFile(df, pathname);
    
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create result set
  outputDataSet <- newInstance(ds, dataFiles);
  setCdf(outputDataSet, cdf);

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})

############################################################################
# HISTORY:
# 2006-12-11
# o Created to immitate the oligo package.
############################################################################
