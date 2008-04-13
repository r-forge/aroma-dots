setMethodS3("readDataFrame", "AffymetrixCdfFile", function(this, units=NULL, fields=NULL, ..., force=FALSE, verbose=FALSE) {
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

  # Argument 'fields':
  if (!is.null(fields)) {
    fields <- Arguments$getCharacters(fields);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading CDF as a data frame");

  key <- list(method="readDataFrame", class=class(this)[1],
              chipType=getChipType(this, fullname=TRUE), 
              units=units, fields=fields, ...);
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


  knownVirtualFields <- c("isPm");
  virtualFields <- intersect(fields, knownVirtualFields);
  verbose && cat(verbose, "Virtual fields: ", paste(virtualFields, collapse=", "));

  verbose2 <- as.integer(isVisible(verbose, -10));
  t0 <- processTime();
  cdfFields <- setdiff(fields, virtualFields);
  if ("isPm" %in% virtualFields) {
    cdfFields <- c(cdfFields, "pbase", "tbase");
  }
  res <- readCdfDataFrame(pathname, units=units, fields=cdfFields, ..., verbose=verbose2);
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


  # Generate virtual fields
  if (length(virtualFields) > 0) {
    if ("isPm" %in% virtualFields) {
      res[["isPm"]] <- with(res, 
        (tbase == "A" & pbase == "T") | 
        (tbase == "T" & pbase == "A") |
        (tbase == "C" & pbase == "G") |
        (tbase == "G" & pbase == "C")
      );
    }
  }

  # Extract fields of interest
  if (!is.null(fields)) {
    res <- res[fields];
  }

  # Save to file cache
  saveCache(res, key=key, dirs=dirs);

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2008-04-13
# o Now readDataFrame() of AffymetrixCdfFile adds "virtual" fields, e.g.
#   the field 'isPm' is inferred and generated from 'pbase' and 'tbase'.
# 2008-04-03
# o Created.
############################################################################
