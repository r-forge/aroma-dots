###########################################################################/**
# @RdocClass DChipQuantileNormalization
#
# @title "The DChipQuantileNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a special @see "QuantileNormalization" that immitates
#  the quantile normalization in the dChip software.
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
setConstructorS3("DChipQuantileNormalization", function(...) {
  extend(QuantileNormalization(...), "DChipQuantileNormalization",
    .exclCells = NULL
  );
})


setMethodS3("getExclCells", "DChipQuantileNormalization", function(this, ..., verbose=FALSE) {
  this$.exclCells;
}, protected=TRUE)


setMethodS3("addExclCells", "DChipQuantileNormalization", function(this, cells, ..., verbose=FALSE) {
  cells <- c(this$.exclCells, cells);
  cells <- unique(cells);
  cells <- sort(cells);
  this$.exclCells <- cells;
}, protected=TRUE)


setMethodS3("excludeChrXFromFit", "DChipQuantileNormalization", function(this, ..., verbose=FALSE) {
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);
  xUnits <- getUnitsOnChromosome(gi, "X");
  cells <- getCellIndices(cdf, units=xUnits);
  cells <- unlist(cells, use.names=FALSE);
  addExclCells(this, cells);
}, protected=TRUE)


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
#   \item{...}{Arguments passed to @see "normalizeQuantileSpline.numeric".}
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
setMethodS3("process", "DChipQuantileNormalization", function(this, ..., force=FALSE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Quantile normalizing (immitating dChip) data set");

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

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Retrieve/calculate the target distribution
  yTarget <- getTargetDistribution(this, verbose=less(verbose));
  yTarget <- sort(yTarget, na.last=TRUE);

  # Get algorithm parameters
  subsetToUpdate <- getSubsetToUpdate(this);

  # Exclude certain cells when *fitting* the normalization function
  excl <- getExclCells(this);
  if (length(excl) > 0) {
    w <- rep(1, length(subsetToUpdate));
    w[excl] <- 0;
    rm(excl);
  } else {
    w <- NULL;
  }

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

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
    x <- getData(df, indices=subsetToUpdate, fields="intensities", verbose=less(verbose,2))$intensities;
    verbose && str(verbose, x);
    verbose && exit(verbose);
  
    x <- normalizeQuantileSpline(x, yTarget, sort=FALSE, w=w, ...);

    # Write normalized data to file
    verbose && enter(verbose, "Writing normalized probe signals");
    # Copy CEL file and update the copy
    verbose && enter(verbose, "Copying source CEL file");
    copyCel(from=getPathname(df), to=pathname, overwrite=!skip);
    verbose && exit(verbose);
    verbose && enter(verbose, "Writing normalized intensities");
    updateCel(pathname, indices=subsetToUpdate, intensities=x);

    rm(x);
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
    verbose && exit(verbose);
  
    # Return new normalized data file object
    dataFiles[[kk]] <- fromFile(df, pathname);
    
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

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
# 2007-03-22
# o Added test code for excluding some cells by giving them weight zero
#   when fitting the normalization function.
# o TO DO: Add code to exclude non-diploid data points from the estimation
#   of the normalization function, e.g. male chromosome X signals should
#   not included if the target distribution was calculated for females only.
#   To simplify it, we could exclude all chromosome X signals regardless
#   of ploidy.  In that we don't have to know the ploidy. The those should
#   also be excluded when estimating the target distribution.
#   Potential problems: If the majority of the signals are from chrX, then
#   it does not work.
# 2006-12-11
# o Created to immitate the oligo package.
############################################################################
