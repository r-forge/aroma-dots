###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod averageQuantile
#
# @title "Gets the average empirical distribution across all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{probes}{An optional @numeric @vector specifying what subset of
#      probes to be used to calculate the empirical distribution. 
#      If @NULL, all probes are used.}
#   \item{...}{Not used.}
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
#   @see "aroma.light::averageQuantile.list"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("averageQuantile", "AffymetrixCelSet", function(this, probes=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'probes':
  probes <- identifyCells(cdf, indices=probes);  # TODO!

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arrays of interest
  arrays <- getNames(this);
  nbrOfChannels <- length(arrays);

  if (is.null(probes)) {
    nbrOfObservations <- nbrOfCells(this);
  } else {
    nbrOfObservations <- length(probes);
  }

  # Construct the sample quantiles
  quantiles <- (0:(nbrOfObservations-1))/(nbrOfObservations-1);

  # Create a vector to hold the target distribution
  xTarget <- vector("double", nbrOfObservations);

  verbose && enter(verbose, "Calculating the average empircal distribution across ", nbrOfChannels, " arrays");

  verbose && printf(verbose, "Number of probes: %d (%.1f%%)\n", 
                   nbrOfObservations, 100*nbrOfObservations/nbrOfCells(cdf));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (cc in 1:nbrOfChannels) {
    verbose && enter(verbose, "Array #", cc);

    verbose && printf(verbose, "reading, ");
    df <- getFile(this, cc);
    Xcc <- getData(df, indices=probes, fields="intensities", ...);
    Xcc <- as.vector(Xcc$intensities);

    # Order and sort the values
    verbose && printf(verbose, "sorting, ");
    Scc <- sort(Xcc);

    # The number of non-NA observations
    nobs <- length(Scc);

    # Has NAs?
    if(nobs < nbrOfObservations) {
      verbose && printf(verbose, "NAs, ");
      tt <- !is.na(Xcc);  # TODO?!? /HB 2006-07-22
      rm(Xcc, tt);

      # Get the sample quantiles for those values
      bins <- (0:(nobs-1))/(nobs-1);

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc'.
      Scc <- approx(x=bins, y=Scc, xout=quantiles, ties="ordered")$y;
      rm(bins);
    } else {
      rm(Xcc);
    }

    # Incremental mean
    verbose && printf(verbose, "summing.");
    xTarget <- xTarget + Scc;
    rm(Scc);

    verbose && exit(verbose);
  }

  xTarget <- xTarget/nbrOfChannels;

  verbose && exit(verbose);

  xTarget;
}) # averageQuantile()




###########################################################################/**
# @RdocMethod normalizeQuantiles
#
# @title "Normalizes samples to have the same empirical distribution"
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
#   \item{subsetToAvg}{The probes to calculate average empirical
#     distribution over.  If a single @numeric in (0,1), then this
#     fraction of all probes will be used.  
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the 
#     average empirical distribution.  
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch 
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{...}{Additional arguments passed to \code{normalizeQuantiles()}.}
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
#   @see "aroma.light::normalizeQuantile.numeric"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("normalizeQuantiles", "AffymetrixCelSet", function(this, outPath=file.path("norm", getChipType(getCdf(this))), xTarget=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, subsetToUpdate=NULL, typesToUpdate=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'outPath':
  outPath <- Arguments$getReadablePathname(outPath, mustExist=FALSE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Get the CDF structure
  cdf <- getCdf(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the average empirical quantiles across all arrays
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  mkdirs(outPath);
  if (is.null(xTarget)) {
    filename <- paste(getChipType(getCdf(this)), "-quantiles.apq", sep="");
    pathname <- filePath(outPath, filename, expandLinks="any");
    verbose && enter(verbose, "Getting average empirical distribution");
    if (isFile(pathname)) {
      verbose && enter(verbose, "Reading saved distribution: ", pathname);
      xTarget <- readApd(pathname)$quantiles;
      verbose && exit(verbose);
    } else {
      probes <- identifyCells(cdf, indices=subsetToAvg, types=typesToAvg,
                                                  verbose=less(verbose));
      verbose && cat(verbose, "Using ", length(probes), " probes");
      verbose && cat(verbose, "Calculating distribution from data");
      xTarget <- averageQuantile(this, probes=probes, 
                                                  verbose=less(verbose));
      verbose && cat(verbose, "Saving distribution: ", pathname);
      writeApd(pathname, data=xTarget, name="quantiles");
    }
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying the probes to be updated");
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate, 
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  verbose && cat(verbose, "Normalizing ", length(subsetToUpdate), " probes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing ", nbrOfArrays(this), " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, "Array #", kk);
    df <- getFile(this, kk);
print(df);
    dataFiles[[kk]] <- normalizeQuantiles(df, outPath=outPath, 
                     xTarget=xTarget, subsetToUpdate=subsetToUpdate, ..., 
                                                  verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  newInstance(this, dataFiles);
})





setMethodS3("transformAffine", "AffymetrixCelSet", function(this, outPath=file.path("transformed", getChipType(getCdf(this))), offset=0, scale=1, subsetToUpdate=NULL, typesToUpdate=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'outPath':
  outPath <- Arguments$getReadablePathname(outPath, mustExist=FALSE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'offset':
  offset <- Arguments$getDouble(offset);

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Identifying the probes to be updated");
  cdf <- getCdf(this);
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate, 
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Transforming ", length(subsetToUpdate),
                               " probes on ", nbrOfArrays(this), " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    df <- getFile(this, kk);
    verbose && enter(verbose, "Array #", kk, " (", getName(df), ")");
    dataFiles[[kk]] <- transformAffine(df, outPath=outPath, 
                  offset=offset, scale=scale, subsetToUpdate=subsetToUpdate,
                                                ..., verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  newInstance(this, dataFiles);
}) # transformAffine()


############################################################################
# HISTORY:
# 2006-09-14
# o Recreated from old AffymetrixDataSet.NORM.R.
# 2006-07-27
# o Added transformAffine().
# o BUG FIX: The 'outPath' argument of normalizeQuantile() in the 
#   AffymetrixDataSet class was not recognized.
# 2006-05-15
# o Extracted from AffymetrixDataSet.R.
# 2006-03-18
# o Added argument 'subset' to calcAvgProbeSignals() & normalizeQuantile().
############################################################################
