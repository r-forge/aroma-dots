setMethodS3("readDataFrame", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  nbrOfUnits <- nbrOfUnits(this);
  if (is.null(units)) {
  } else {
    # Validate unit indices
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits));
    nbrOfUnits <- length(units);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading CDF as a data frame");

  key <- list(method="readDataFrame", class=class(this)[1],
              chipType=getChipType(this, fullname=TRUE), units=units, ...);
  dirs <- c("aroma.affymetrix", getChipType(this, fullname=TRUE));
  if (!force) {
    res <- loadCache(key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Found results cached on file");
      verbose && exit(verbose);
      return(res);
    }
  }

  if (!exists("readCdfDataFrame", mode="function")) {
    throw("This method requires readCdfDataFrame() in affxparser.");
  }

  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);


  verbose2 <- as.integer(isVisible(verbose, -10));
  t0 <- processTime();
  res <- readCdfDataFrame(pathname, units=units, ..., verbose=verbose2);
  t1 <- processTime();

  verbose && cat(verbose, "Read data:");
  verbose && str(verbose, res);

  # Make nucleotide bases in upper case.
  for (field in c("pbase", "tbase")) {
    res[[field]] <- toupper(res[[field]]);
  }

  if (verbose) {
    dt <- (t1-t0)[3];
    cat(verbose, "Total reading/processing time:");
    printf(verbose, "Total time: %.0f secs = %.2f mins = %.2f hours\n", 
                                                        dt, dt/60, dt/3600);
    printf(verbose, "Time/unit: %.2f ms = %.2f secs\n", 
                                         1000*dt/nbrOfUnits, dt/nbrOfUnits);
  }

  # Save to file cache
  saveCache(res, key=key, dirs=dirs);

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2008-04-03
# o Created.
############################################################################
