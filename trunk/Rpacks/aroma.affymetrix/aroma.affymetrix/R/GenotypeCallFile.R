###########################################################################/**
# @RdocClass GenotypeCallFile
#
# @title "The GenotypeCallFile class"
#
# \description{
#  @classhierarchy
#
#  An GenotypeCallFile object represents parameter estimates.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# @visibility "private"
#*/###########################################################################
setConstructorS3("GenotypeCallFile", function(..., cdf=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cdf)) {
    if (is.character(cdf)) {
      cdf <- AffymetrixCdfFile$fromChipType(cdf);
    }

    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  extend(AffymetrixFile(...), "GenotypeCallFile",
    cdf = cdf
  )
})

setMethodS3("create", "GenotypeCallFile", function(static, filename, path="calls", cdf, overwrite=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(cdf)) {
    cdf <- AffymetrixCdfFile$fromChipType(cdf);
  } else if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }
  nbrOfUnits <- nbrOfUnits(cdf);
  nbrOfUnits <- Arguments$getInteger(nbrOfUnits, range=c(1,Inf));

  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);

  calls <- FileByteVector(pathname, length=nbrOfUnits);
  on.exit(close(calls));

  newInstance(static, pathname, cdf=cdf);
}, static=TRUE, protected=TRUE)



setMethodS3("getCdf", "GenotypeCallFile", function(this, ...) {
  this$cdf;
})

setMethodS3("[", "GenotypeCallFile", function(this, i, drop=TRUE, ...) {
  res <- readUnits(this, units=i, ...);
  if (drop && length(res) == 1)
    res <- drop(res);
  res;
})

setMethodS3("[<-", "GenotypeCallFile", function(this, i, value) {
  writeUnits(this, units=i, calls=value);
})

setMethodS3("readUnits", "GenotypeCallFile", function(this, units=NULL, ...) {
  # Argument 'units':
  cdf <- getCdf(this);
  if (is.null(units)) {
    units <- 1:nbrOfUnits(cdf);
  } else if (is.character(units)) {
    units <- indexOf(cdf, units);
  }

  # Open data file
  pathname <- getPathname(this);
  fdata <- FileByteVector(pathname);
  on.exit(close(fdata));

  # Retrieve data
  calls <- factor(fdata[units], levels=0:4);
  levels(calls) <- c("NotCalled", "AA", "AB", "BB", "NoCall");

  calls;
})


setMethodS3("updateUnits", "GenotypeCallFile", function(this, units=NULL, calls, ...) {
  # Argument 'units':
  cdf <- getCdf(this);
  if (is.null(units)) {
    units <- 1:nbrOfUnits(cdf);
  } else if (is.character(units)) {
    units <- indexOf(cdf, units);
  }

  # Argument 'calls':
  if (is.numeric(calls)) {
    calls <- as.integer(calls);
  } else if (is.factor(calls)) {
    calls <- as.integer(calls);
  } else if (is.character(calls)) {
    calls <- match(calls, c("AA", "AB", "BB"));
    calls[is.na(calls)] <- 0;
  }
  if (length(calls) != length(units)) {
    throw("Argument 'units' and 'calls' are of different lengths: ", 
                                     length(calls), " != ", length(units));
  }

  # Open data file
  pathname <- getPathname(this);
  fdata <- FileByteVector(pathname);
  on.exit(close(fdata));

  # Update data
  fdata[units] <- calls;

  invisible(this);
})


setMethodS3("importFromCrlmmFile", "GenotypeCallFile", function(this, filename="calls.txt", path=NULL, sample=getName(this), ..., verbose=FALSE) {
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing genotype calls from CRLMM output file.");
  verbose && cat(verbose, "Pathname: ", pathname);

  # Read column header
  fh <- file(pathname, open="r");
  on.exit(close(fh));

  # Get the sample names
  header <- readLines(con=fh, n=1);
  sep <- " ";
  if (regexpr("\t", header) != -1)
    sep <- "\t";
  header <- unlist(strsplit(header, split=sep));
  header <- trim(header);
  header <- gsub("^\"", "", header);
  header <- gsub("\"$", "", header);
  header <- gsub("[.](c|C)(e|E)(l|L)$", "", header);
  nbrOfSamples <- length(header);

  verbose && printf(verbose, "Detected %d samples.", nbrOfSamples);

  if (is.character(sample)) {
    sampleIdx <- grep(sample, header);
    if (length(sampleIdx) == 0) {
      throw("Could not import genotype calls. No sample '", sample, 
                                          "' in call file: ", pathname);
    } else if (length(sampleIdx) > 2) {
      throw("Could not import genotype calls. More than one match for sample '", sample, "' in call file: ", pathname);
    }
    sample <- sampleIdx;
  }
  verbose && printf(verbose, "Requested sample in column %d.", sample);
  sample <- Arguments$getIndex(sample, range=c(1, nbrOfSamples));

  verbose && enter(verbose, "Reading data");
  colClasses <- c("character", rep("NULL", nbrOfSamples));
  colClasses[sample+1] <- "integer";
  calls <- read.table(fh, colClasses=colClasses, row.names=1, sep=sep, header=FALSE);
  ncalls <- nrow(calls);
  verbose && cat(verbose, "Number of rows read: ", ncalls);
  verbose && exit(verbose);

  # Sort calls according to CDF
  verbose && enter(verbose, "Sort according to CDF file");
  cdf <- getCdf(this);
  verbose && cat(verbose, "Number of units in CDF: ", nbrOfUnits(cdf));
  units <- indexOf(cdf, names=rownames(calls));
  nbad <- which(is.na(units));
  
  # Check for units not in CDF and exclude them
  if (length(nbad) == 0) {
    verbose && cat(verbose, "All ", length(units), " SNPs read are available in the CDF file.");
  } else {
    warning(sprintf("%d (%%.1f) of the SNPs in the call file was not found in the CDF file: %s", nbad, 100*nbad/ncalls, pathname));
    units <- units[-nbad];
  }
  verbose && exit(verbose);

  calls <- calls[,1];

  # Update
  updateUnits(this, units=units, calls=calls, verbose=less(verbose));

  verbose && exit(verbose);

  invisible(length(units));
})



############################################################################
# HISTORY:
# 2006-10-01
# o Can now import a single CRLMM column.
# 2006-09-30
# o Created.
############################################################################
