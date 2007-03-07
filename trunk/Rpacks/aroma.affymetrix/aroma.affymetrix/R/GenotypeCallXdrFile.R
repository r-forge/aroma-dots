###########################################################################/**
# @RdocClass GenotypeCallXdrFile
#
# @title "The GenotypeCallXdrFile class"
#
# \description{
#  @classhierarchy
#
#  The abstract GenotypeCallXdrFile class represents a file containing genotype
#  calls for a given genotype parameter estimates.
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
setConstructorS3("GenotypeCallXdrFile", function(...) {
  this <- extend(GenotypeCallFile(...), "GenotypeCallXdrFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("getCdf", "GenotypeCallXdrFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    path <- getPath(this);
    chipType <- basename(path);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("fromFile", "GenotypeCallXdrFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  df <- newInstance(static, pathname, ...);

  df;
}, static=TRUE);


setMethodS3("readUnits", "GenotypeCallXdrFile", function(this, units=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  units <- convertUnits(cdf, units=units, keepNULL=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open data file
  pathname <- getPathname(this);
  calls <- loadObject(pathname);
  if (!is.null(units)) {
    calls <- calls[units];
  }

  calls;
})


setMethodS3("updateUnits", "GenotypeCallXdrFile", function(this, units=NULL, calls, ...) {
  throw("Not implemented.");
})




############################################################################
# HISTORY:
# 2006-12-13
# o Revived.
# 2006-10-01
# o Can now import a single CRLMM column.
# 2006-09-30
# o Created.
############################################################################
