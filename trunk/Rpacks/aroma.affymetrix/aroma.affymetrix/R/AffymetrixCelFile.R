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
setMethodS3("readUnits", "AffymetrixCelFile", function(this, units=NULL, cdf=NULL, ..., stratifyBy=NULL) {
  if (is.null(cdf)) {
    suppressWarnings({
      cdf <- readUnits(getCdf(this), units=units, stratifyBy=stratifyBy);
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



setMethodS3("[", "AffymetrixCelFile", function(this, units=NULL, drop=FALSE) {
  data <- readUnits(this, units=units);
  if (drop && length(data) == 1)
    data <- data[[1]];
  data;
})

setMethodS3("[[", "AffymetrixCelFile", function(this, unit=NULL) {
  this[units=unit, drop=TRUE];
})


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


setMethodS3("getRectangle", "AffymetrixCelFile", function(this, xrange=c(0,Inf), yrange=c(0,Inf), fields=c("intensities", "stdvs", "pixels"), ...) {
  readCelRectangle(this$.pathname, xrange=xrange, yrange=yrange, readIntensities=("intensities" %in% fields), readStdvs=("stdvs" %in% fields), readPixels=("pixels" %in% fields));
})



###########################################################################/**
# @RdocMethod image270
#
# @title "Displays all or a subset of the data spatially"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{xrange}{A @numeric @vector of length two giving the left and right
#          coordinates of the cells to be returned.}
#   \item{yrange}{A @numeric @vector of length two giving the top and bottom
#          coordinates of the cells to be returned.}
#   \item{...}{Additional arguments passed to @seemethod "readRectangle",
#      but also @see "graphics::image".}
#   \item{field}{The data field to be displayed.}
#   \item{col}{The color map to be used.}
#   \item{main}{The main title of the plot.}
# }
#
# \value{
#  Returns the (270-degrees rotated) data @matrix.
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
setMethodS3("image270", "AffymetrixCelFile", function(this, xrange=c(0,Inf), yrange=c(0,Inf), ..., field=c("intensities", "stdvs", "pixels"), col=gray.colors(256), main=getName(this)) {
  rotate270 <- function(x, ...) {
    x <- t(x)
    nc <- ncol(x)
    if (nc < 2) return(x)
    x[,nc:1,drop=FALSE]
  }

  # Argument 'field':
  field <- match.arg(field);

  suppressWarnings({
    cel <- getRectangle(this, xrange=xrange, yrange=yrange, fields=field, ...);
  })

  # Display the rectangle
  y <- cel[[field]];

  suppressWarnings({
    image(rotate270(y), col=col, ..., axes=FALSE, main=main);
  })

  if (is.null(xrange)) 
    xrange <- c(0,ncol(y)-1);
  if (is.null(yrange)) 
    yrange <- c(0,nrow(y)-1);

  cdf <- getCdf(this);
  dim <- paste(getDimension(cdf), collapse="x");
  label <- sprintf("Chip type: %s [%s]", getChipType(cdf), dim);
  text(x=0, y=0, labels=label, adj=c(0,1.2), cex=0.8, xpd=TRUE)
  label <- sprintf("(%d,%d)", as.integer(xrange[1]), as.integer(yrange[1]));
  text(x=0, y=1, labels=label, adj=c(0,-0.7), cex=0.8, xpd=TRUE)
  label <- sprintf("(%d,%d)", as.integer(xrange[2]), as.integer(yrange[2]));
  text(x=1, y=0, labels=label, adj=c(1,1.2), cex=0.8, xpd=TRUE)

  # Return the plotted data.
  invisible(y);
})




###########################################################################/**
# @RdocMethod writeSpatial
#
# @title "Writes a spatial image of the CEL data to an image file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the image file.}
#   \item{path}{The path to the image file.}
#   \item{...}{Additional arguments passed to @seemethod "readRectangle".}
#   \item{field}{The data field to be displayed.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns (invisibly) the pathname of the generated image file.
# }
#
# \section{Details}{
#   If any other image format than the pixmap *.pgm format is requested,
#   an external image coverted is used to convert from the PGM format.
#   You may specify the absolute pathname to such an image converted by:
#   \code{options(imageConverter="C:/Program Files/ImageMagick/bin/convert")}.
# }
#
# @author
#
# \seealso{
#   @seemethod "image270".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("writeSpatial", "AffymetrixCelFile", function(this, filename=sprintf(fmtstr, getName(this)), path=filePath("figures", getChipType(getCdf(this))), fmtstr="%s-spatial.png", ..., field=c("intensities", "stdvs", "pixels"), verbose=FALSE) {
  require(R.image) || throw("Package R.image not loaded.");

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  verbose && enter(verbose, "Read CEL data");
  suppressWarnings({
    cel <- getRectangle(this, fields=field, ...);
  })
  verbose && exit(verbose);

  verbose && enter(verbose, "Create in memory image");
  # Get the image
  img <- cel[[field]];
  
  # Create 256-levels gray-scale image
  img <- sqrt(img) - 1;
  img <- GrayImage(img);
  verbose && exit(verbose);

  verbose && enter(verbose, "Write image to temporary PGM file");
  pgmfile <- paste(pathname, ".pgm", sep="");
  on.exit({
    # Remove the PGM file if the final file was successfully create.
    if (isFile(pathname))
      file.remove(pgmfile)
  });
  write(img, pgmfile);
  verbose && exit(verbose);

  # Done?
  if (regexpr("[.]pgm$", pathname) != -1)
    return(invisible(pathname));

  # Convert PGM to final image format using ImageMagick
  verbose && enter(verbose, "Convert to final file using ImageMagick");
  converter <- getOption("imageConverter");
  if (is.character(converter)) {
    converter <- function(srcfile, destfile, format, options=NULL, ...) {
      res <- system(paste(converter, options, pgmfile, pathname));
      (res == 0);
    }
  }

  suppressWarnings({
    converter(pgmfile, pathname, ...);
  })

  verbose && exit(verbose);

  invisible(pathname);
})


############################################################################
# HISTORY:
# 2006-08-25
# o Added image270() and writeSpatial().
# o Added methods "[" and "[[" mapped to readUnits().
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
