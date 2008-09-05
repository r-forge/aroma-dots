# SPEEDUP TRICK: Update CN units (which are all single probe units) by copying
setMethodS3("fitCnProbes", "UnitModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constants
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .Machine$float.eps <- sqrt(.Machine$double.eps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Estimating single-probe CN units");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  params <- getParameters(this);
  shift <- params$shift;

  verbose && enter(verbose, "Getting chip-effect set");
  # Get chip-effect set
  ces <- getChipEffectSet(this, verbose=verbose);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying CN units");
  cdf <- getCdf(ds);
  units <- whichVector(getUnitTypes(cdf) == 5);
  verbose && str(verbose, units);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying cell indices for CN units");
  cells <- getCellIndices(cdf, units=units, useNames=FALSE, unlist=TRUE);
  verbose && str(verbose, cells);
  # Sanity check
  if (length(cells) != length(units)) {
    throw("Detected units with more than one cell: ", getChipType(cdf));
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying cell indices for estimates");
  cdfM <- getCdf(ces);
  cellsM <- getCellIndices(cdfM, units=units, useNames=FALSE, unlist=TRUE);
  verbose && str(verbose, cellsM);
  # Sanity check
  if (length(cellsM) != length(units))
    throw("Detected units with more than one cell: ", getChipType(cdfM));
  if (length(cellsM) != length(cells)) {
    throw("Input 'cells' and output 'cellsM' are of different lengths: ", 
                                  length(cellsM), " != ", length(cells));
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # "Fitting"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting ", length(ds), " arrays");
  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), length(ds)));

    verbose && enter(verbose, "Reading signals");
    cef <- getFile(ces, kk);
    y <- extractMatrix(df, cells=cells, drop=TRUE);
    stopifnot(length(y) == length(cells));
    verbose && str(verbose, y);
    verbose && exit(verbose);

    # Shift?
    if (shift != 0) {
      verbose && enter(verbose, "Shifting signals");
      y <- y + shift;
      verbose && str(verbose, y);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Transforming signals to estimates");
    sdTheta <- .Machine$float.eps;  # Smallest float > 0.
    data <- data.frame(cell=cellsM, theta=y, sdTheta=sdTheta, outliers=FALSE);
    rm(y);
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing estimates");
    updateDataFlat(cef, data=data, verbose=log);
    rm(data);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(units);
}, private=TRUE) # fitCnProbes()


############################################################################
# HISTORY:
# 2008-09-05
# o Created.
############################################################################
