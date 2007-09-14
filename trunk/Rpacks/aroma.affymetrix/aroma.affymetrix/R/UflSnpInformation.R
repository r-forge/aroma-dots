setConstructorS3("UflSnpInformation", function(..., .ufl=NULL, .verify=TRUE) {
  this <- extend(SnpInformation(...), "UflSnpInformation",
    .ufl = .ufl
  );
  if (.verify) {
    if (!is.null(getPathname(this)))
      verify(this);
  }
  this;
})

setMethodS3("getAromaUflFile", "UflSnpInformation", function(this, ..., force=FALSE) {
  ufl <- this$.ufl;
  if (force || is.null(ufl)) {
    ufl <- AromaUflFile(getPathname(this), ...);
    this$.ufl <- ufl;
  }
  ufl;
}, protected=TRUE);


setMethodS3("findByChipType", "UflSnpInformation", function(static, ...) {
  AromaUflFile$findByChipType(...);
}, static=TRUE, protected=TRUE)


setMethodS3("fromChipType", "UflSnpInformation", function(static, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ufl <- AromaUflFile$fromChipType(...);
  pathname <- getPathname(ufl);

  verbose && enter(verbose, "Instantiating ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);

  res <- newInstance(static, filename=pathname, path=NULL, .ufl=ufl, .verify=FALSE);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})


setMethodS3("verify", "UflSnpInformation", function(this, ...) {
  tryCatch({
    df <- readData(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the UFL SNP information file (", 
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readData", "UflSnpInformation", function(this, nrow=NULL, ..., verbose=FALSE) {
  verbose && enter(verbose, "Reading data from UFL file");

  ufl <- getAromaUflFile(this);
  verbose && print(verbose, ufl, level=-20);

  if (is.null(nrow)) {
    verbose && cat(verbose, "Reading all ", nbrOfUnits(ufl), " units");
    res <- ufl[,,drop=FALSE];
  } else {
    units <- 1:nrow;
    verbose && cat(verbose, "Reading ", length(units), " units");
    res <- ufl[units,,drop=FALSE];
  }

  colnames(res) <- c("fragmentLength");
  verbose && exit(verbose);

  res;
})

setMethodS3("getDataColumns", "UflSnpInformation", function(this, ...) {
  ufl <- getAromaUflFile(this);
  nbrOfColumns <- nbrOfColumns(ufl);
  names <- c("fragmentLength", rep(c("fragmentLength"), nbrOfColumns-1));
  if (nbrOfColumns > 1)
    names[-1] <- sprintf("%s,%02d", names[-1], 2:nbrOfColumns);
  names;
}, private=TRUE)

setMethodS3("getFields", "UflSnpInformation", function(this, ...) {
  getDataColumns(this, ...);
})

setMethodS3("nbrOfUnits", "UflSnpInformation", function(this, ...) {
  ufl <- getAromaUflFile(this);
  nbrOfUnits(ufl);
})


setMethodS3("getData", "UflSnpInformation", function(this, units=NULL, fields=getDataColumns(this), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving SNP information from file");

    # Now read the genome information data
    ufl <- getAromaUflFile(this);
    cc <- match(fields, getDataColumns(this));
    missing <- fields[is.na(cc)];
    if (length(missing)) {
      throw("Unknown fields: ", paste(missing, collapse=", "));
    }
  
    verbose && enter(verbose, "Reading SNP information data");
    data <- ufl[,,drop=FALSE];
    colnames(data) <- getDataColumns(this);
    verbose && str(verbose, data);
    verbose && exit(verbose);

    # Store in cache
    this$.data <- data;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    
    verbose && exit(verbose);
  }

  # Subset by unit?
  if (!is.null(units)) {
    # Map the unit indicies to the row names
    data <- data[units,,drop=FALSE];
  }

  # Stratify by field values?
  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      # Get the values to be stratified upon.
      values <- data[,key,drop=FALSE];

      # Get the test (value or function)
      test <- args[[key]];
      test <- na.omit(test);
      if (is.function(test)) {
        keep <- test(values);
      } else {
        keep <- (values == test);
        keep <- (keep & !is.na(keep));
      }
      data <- data[keep,,drop=FALSE];
    }
    rm(keep);
  }

  # Reorder?
  if (!is.null(orderBy)) {
    o <- do.call("order", args=as.list(data[,orderBy,drop=FALSE]));
    data <- data[o,,drop=FALSE];
    rm(o);
  }

  # Extract a subset of fields?
  if (!is.null(fields))
    data <- data[,fields, drop=FALSE];

  data;
})


setMethodS3("getFragmentStarts", "UflSnpInformation", function(this, ...) {
  throw("Not supported.");
})


setMethodS3("getFragmentStops", "UflSnpInformation", function(this, ...) {
  throw("Not supported.");
})


############################################################################
# HISTORY:
# 2007-09-11
# o Created from DChipSnpInformation.R.
############################################################################  
