setMethodS3("getCalculateResidualsFunction", "ProbeLevelModel", function(static, ...) {
  function(y, yhat) {
    y-yhat;
  }
}, static=TRUE, protected=TRUE)


setMethodS3("calculateResiduals", "ProbeLevelModel", function(this, units=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ces <- getChipEffects(this);
  paf <- getProbeAffinities(this);
  rs <- getResidualSet(this);
  nbrOfArrays <- nbrOfArrays(ces);

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  verbose && enter(verbose, "Calculating PLM residuals");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data and parameter objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  if (is.null(ds)) {
    throw("No data set specified for PLM: ", getFullName(this));
  }

  # Get the function how to calculate residuals.
  # Default is eps = y - yhat, but for instance RMA uses eps = y/yhat.
  calculateEps <- getCalculateResidualsFunction(this);

  cdf <- getCdf(ds);
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  cdfData <- NULL;
  chipType <- getChipType(cdf);
  key <- list(method="calculateResiduals", class=class(this)[1], 
              chipType=chipType, params=getParameters(this),
              units=units);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    cdfData <- loadCache(key, dirs=dirs);
    if (!is.null(cdfData))
      verbose && cat(verbose, "Found indices cached on file");
  }

  if (is.null(cdfData)) {
    verbose && enter(verbose, "Retrieving CDF cell indices");
    cdfUnits <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculate group sizes");
    unitGroupSizes <- applyCdfGroups(cdfUnits, lapply, FUN=function(group) {
      length(.subset2(group, 1));
    });
    unitGroupSizes <- unlist(unitGroupSizes, use.names=FALSE);
    verbose && exit(verbose);

    cells <- unlist(cdfUnits, use.names=FALSE);

    verbose && enter(verbose, "Retrieving CDF cell indices for chip effects");
    cdfUnits <- getCellIndices(ces, units=units, verbose=less(verbose));
    cells2 <- unlist(cdfUnits, use.names=FALSE);
    verbose && exit(verbose);

    rm(cdfUnits); # Not needed anymore

    cdfData <- list(unitGroupSizes=unitGroupSizes, cells=cells, cells2=cells2);
    verbose && enter(verbose, "Saving to file cache");
    saveCache(cdfData, key=key, dirs=dirs);
    verbose && exit(verbose);
  } else {
    unitGroupSizes <- cdfData$unitGroupSizes;
    cells <- cdfData$cells;
    cells2 <- cdfData$cells2;
    rm(cdfData);
  }

  # Optimized reading order
  o <- .Internal(qsort(cells, TRUE));
  cells <- o$x;
  o <- o$ix;
  oinv <- .Internal(qsort(o, TRUE))$ix;
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Retrieving probe-affinity estimates");
  phi <- getData(paf, indices=cells, fields="intensities")$intensities[oinv];
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);

  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    cef <- getFile(ces, kk);
    rf <- getFile(rs, kk);

    verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, getName(df)));

    filename <- sprintf("%s,residuals.cel", getFullName(df));
    pathname <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "Pathname: ", pathname);
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already calculated.");
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Retrieving probe intensity data");
    y <- getData(df, indices=cells, fields="intensities")$intensities[oinv];
    verbose && exit(verbose);

    verbose && enter(verbose, "Retrieving chip-effect estimates");
    theta <- getData(cef, indices=cells2, fields="intensities")$intensities;
    theta <- rep(theta, times=unitGroupSizes);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating residuals");
    yhat <- phi * theta;
    eps <- calculateEps(y, yhat);  # Model class specific.
    verbose && str(verbose, eps);
    verbose && exit(verbose);
    rm(y, yhat, theta);

    verbose && enter(verbose, "Storing residuals");
    tryCatch({
      # Copy CEL file and update the copy
      verbose && enter(verbose, "Copying source CEL file");
      copyCel(from=getPathname(df), to=pathname, overwrite=force);
      verbose && exit(verbose);
      verbose && enter(verbose, "Writing normalized intensities");
      updateCel(getPathname(rf), indices=cells, intensities=eps[o]);
      verbose && exit(verbose);
    }, interrupt = function(intr) {
      verbose && print(verbose, intr);
      file.remove(pathname);
    }, error = function(ex) {
      verbose && print(verbose, ex);
      file.remove(pathname);
    })
    verbose && exit(verbose);

#    verbose && enter(verbose, "Verifying");
#    eps2 <- getData(rf, indices=cells, fields="intensities")$intensities[oinv];
#    stopifnot(all.equal(eps, eps2, tolerance=.Machine$double.eps^0.25));
#    verbose && exit(verbose);
#    rm(eps2);

    rm(eps);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (kk ...)
  rm(cells, phi, unitGroupSizes);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  invisible(rs);
})


##########################################################################
# HISTORY:
# 2007-02-14 HB + KS
# o Now residuals can be calculated differently for different PLM classes.
#   This is done by overriding static getCalculateResidualsFunction().
# 2007-02-12 HB
# o Rewritten from KS:s code.
##########################################################################
