###########################################################################/**
# @RdocClass AffymetrixCelFile
#
# @title "The AffymetrixCelFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCelFile object represents a single Affymetrix CEL file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
#   \item{cdf}{An optional @see "AffymetrixCdfFile" making it possible to
#     override the default CDF file as specified by the CEL file header.  
#     The requirement is that its number of cells must match that of
#     the CEL file.
#     If @NULL, the CDF structure is inferred from the the chip type
#     as specified in the CEL file header.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \section{For developers}{
#   If you subclass this class, please make sure to query the
#   @see "AffymetrixCdfFile" object (see @seemethod "getCdf") whenever
#   querying CDF information.  Do not use the CDF file inferred from the
#   chip type in CEL header, unless you really want it to be hardwired that
#   way, otherwise you will break to possibility to override the CDF 
#   structure.
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixCelFile", function(..., cdf=NULL) {
  this <- extend(AffymetrixFile(...), "AffymetrixCelFile",
    "cached:.header" = NULL,
    .cdf = NULL
  )

  if (!is.null(cdf))
    setCdf(this, cdf);

  this;
})


###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixCelFile object from a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelFile" object.  
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AffymetrixCelFile", function(static, filename, path=NULL, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  # WORKAROUND: Currently the affxparser code crash R if the file is not
  # a valid CEL file.  The best we can do now is to test against the
  # filename.
  isCel <- (regexpr("[.](c|C)(e|E)(l|L)$", pathname) != -1);
  if (!isCel) {
    throw("Could not read CEL file. File format error: ", pathname);
  }

  # Try to read the header assuming
  header <- readCelHeader(pathname);

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE)




###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "AffymetrixCelFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chiptype <- getHeader(this)$chiptype;
    cdf <- this$.cdf <- AffymetrixCdfFile$fromChipType(chiptype);
  }
  cdf;
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "AffymetrixCelFile", function(this, cdf, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile object: ", 
                                                               class(cdf)[1]);
  }

  # Assure that the CDF is compatible with the CEL file
  if (nbrOfCells(cdf) != getHeader(this)$total) {
    throw("The specified CDF structure is not compatible with the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", getHeader(this)$total);
  }

  this$.cdf <- cdf;
  invisible(this);
})



###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list structure as returned by @see "affxparser::readCelHeader".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AffymetrixCelFile", function(this, ...) {
  header <- this$.header;
  if (is.null(header))
    header <- this$.header <- readCelHeader(this$.pathname);
  header;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type for this CEL file"
#
# \description{
#  @get "title" \emph{according} to the @see "AffymetrixCdfFile" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getChipType", "AffymetrixCelFile", function(this, ...) {
  getChipType(getCdf(this));
}, protected=TRUE)



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads CEL data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{cdf}{An alternative CDF structure to be used.  This overrides
#     the \code{units} arguments.}
#   \item{...}{Arguments passed to \code{getUnits()} of the 
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the @see "affxparser::readCelUnits" methods.}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::readCelUnits" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "AffymetrixCelFile", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf)) {
    suppressWarnings({
      cdf <- getUnits(getCdf(this), units=units, ...);
    });
  }
  suppressWarnings({
    readCelUnits(this$.pathname, cdf=cdf, ...);
  })
}, protected=TRUE);


###########################################################################/**
# @RdocMethod updateUnits
#
# @title "Updates CEL data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units.
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @list structure consisting of grouped data similar to 
#      what @seemethod "readUnits" returns.}
#   \item{...}{Additional arguments passed to the
#      @see "affxparser::updateCelUnits" methods.}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::updateCelUnits" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("updateUnits", "AffymetrixCelFile", function(this, data, ...) {
  updateCelUnits(this$.pathname, data=data, ...);
}, protected=TRUE);




setMethodS3("getIntensities", "AffymetrixCelFile", function(this, indices=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  nbrOfCells <- nbrOfCells(getCdf(this));
  if (is.null(indices)) { 
  } else {
    indices <- Arguments$getIndices(indices, range=c(1,nbrOfCells));
    nbrOfCells <- length(indices);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  pathname <- this$.pathname;
  cVerbose <- -(as.numeric(verbose) + 1);
  y <- readCel(pathname, indices=indices, 
          readHeader=FALSE, 
          readIntensities=TRUE, readStdvs=FALSE, 
          readPixels=FALSE, readXY=FALSE, 
          readOutliers=FALSE, readMasked=FALSE,
          ...,
          verbose=cVerbose)$intensities;
   y;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2006-08-24
# o Added the option to specify an 'cdf' object, making it possible to 
#   override the default CDF file according to the CEL header.  It is
#   important that all methods queries the AffymetrixCdfFile object
#   from getCdf() and not the one through the CEL header!
# o Added most Rdoc comments.
# 2006-07-21
# o Added readUnits().
# 2006-07-05
# o BUG FIX/WORKAROUND: Currently the affxparser code crash R if the file 
#   is not a valid CEL file.  The best we can do now is to test that the
#   filename has suffix *.CEL.
# 2006-05-30
# o Added fromFile().
# 2006-03-30
# o Updated according to affxparser.
# 2006-03-23
# o Moved all SNP related methods into the new class AffymetrixSnpCelFile.
# 2006-03-18
# o Made probe indices one-based.
# 2006-03-04
# o Added support for remapping in readIntensities().  This is currently
#   not used for CEL files (only APD files), but was added for the future.
# 2006-03-02
# o Created.
############################################################################
