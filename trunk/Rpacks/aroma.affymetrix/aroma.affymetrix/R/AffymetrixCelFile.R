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
#
# \seealso{
#   An object of this class is typically part of an @see "AffymetrixCelSet".
# }
#*/###########################################################################
setConstructorS3("AffymetrixCelFile", function(..., cdf=NULL) {
  this <- extend(AffymetrixFile(...), "AffymetrixCelFile",
    "cached:.header" = NULL,
    "cached:lastPlotData" = NULL,
    .cdf = NULL
  )

  if (!is.null(cdf))
    setCdf(this, cdf);

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    parseTagsAsAttributes(this);

  this;
})

setMethodS3("clearCache", "AffymetrixCelFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".header", ".lastPlotData")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("clone", "AffymetrixCelFile", function(this, ..., verbose=TRUE) {
  # Clone itself (and clear the cached fields)
  object <- NextMethod("clone", clear=TRUE, ...);

  # Clone the CDF here.
  if (!is.null(object$.cdf))
    object$.cdf <- clone(object$.cdf);

  object;
})


setMethodS3("as.character", "AffymetrixCelFile", function(this, ...) {
  s <- NextMethod("as.character", ...);
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("Timestamp: %s", as.character(getTimestamp(this))));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getIdentifier", "AffymetrixCelFile", function(this, ..., force=FALSE) {
  identifier <- this$.identifier;
  if (force || is.null(identifier)) {
    # Get header
    hdr <- getHeader(this);
    # Get subset of data
    nbrOfCells <- hdr$total;
    mid <- nbrOfCells %/% 2;
    subset <- seq(from=mid - 500, to=mid + 500);
    data <- getData(this, indices=subset);
    key <- list(hdr=hdr, data=data);
    id <- digest(key);
    this$.identifier <- id;
  }
  identifier;
}, private=TRUE)




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
setMethodS3("fromFile", "AffymetrixCelFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
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


  # WORKAROUND: Currently the affxparser code crash R if the file is not
  # a valid CEL file.  The best we can do now is to test against the
  # filename.
  isCel <- (regexpr("[.](c|C)(e|E)(l|L)$", pathname) != -1);
  if (!isCel) {
    throw("Could not read CEL file. Filename format error: ", pathname);
  }

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("createFrom", "AffymetrixCelFile", function(this, filename, path=NULL, clear=TRUE, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Don't create if already exists
  if (isFile(pathname)) {
    res <- newInstance(this, pathname);
    setCdf(res, getCdf(this));
    return(res);
  }

  originalCdf <- getCdf(this);
  
  res <- copyFile(this, filename=pathname, path=NULL, verbose=less(verbose));

  if (clear) {
    clearData(res, ..., .forSure=TRUE, verbose=less(verbose));
  }

  setCdf(res, originalCdf);
  res;
}, private=TRUE)



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
    chipType <- getHeader(this)$chiptype;
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
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
setMethodS3("setCdf", "AffymetrixCelFile", function(this, cdf, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'cdf':
    if (!inherits(cdf, "AffymetrixCdfFile")) {
      throw("Argument 'cdf' is not an AffymetrixCdfFile object: ", 
                                                                 class(cdf)[1]);
    }
  
    # Assure that the CDF is compatible with the CEL file
    if (nbrOfCells(cdf) != nbrOfCells(this)) {
      throw("The specified CDF structure is not compatible with the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(this));
    }

    # Nothing to do?
#    oldCdf <- getCdf(this);
#    if (equals(cdf, oldCdf))
#      return(invisible(this));
  }

  # Have to clear the cache 
  clearCache(this);

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
}, private=TRUE)


setMethodS3("getHeaderV3", "AffymetrixCelFile", function(this, ...) {
  # Get the CEL header
  header <- getHeader(this);

  # Get the CEL v3 header
  header <- header$header;

  # Parse it
  header <- unlist(strsplit(header, split="\n"));
  header <- trim(header);
  header <- header[nchar(header) > 0];
  header <- strsplit(header, split="=");
  names <- sapply(header, FUN=function(s) s[1]);
  header <- lapply(header, FUN=function(s) s[-1]);
  names(header) <- names;

  header;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getTimestamp
#
# @title "Gets the timestamp in the CEL header"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{format}{The default format string for parsing the time stamp.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a \code{POSIXct} object.
#  The parsed string containing the timestamp is returned as 
#  attribute \code{text}.
# }
#
# @author
#
# \seealso{
#   Internally, @see "base::strptime" is used to parse the time stamp.
#   @see "base::DateTimeClasses" for more information on \code{POSIXct}.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getTimestamp", "AffymetrixCelFile", function(this, format="%m/%d/%y %H:%M:%S", ...) {
  # Argument 'format':
  format <- Arguments$getCharacter(format);


  # Get the CEL v3 header of the CEL header
  header <- getHeaderV3(this);

  # Get the DAT header
  header <- header$DatHeader;

  # Find the element with a date. It is part of the same string as the
  # one containing the chip type.  Get the chip type from the header.
  chipType <- getHeader(this)$chiptype;
  pattern <- sprintf(" %s.1sq ", chipType);
  header <- grep(pattern, header, value=TRUE);

  # Extract the date timestamp
  pattern <- ".*([01][0-9]/[0-3][0-9]/[0-9][0-9] [0-2][0-9]:[0-5][0-9]:[0-5][0-9]).*";
  timestamp <- gsub(pattern, "\\1", header);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alternative:
  # Could use a pattern, but if a different timestamp than the American is 
  # used, this wont work.  Instead assume a fixed location.
  # From the DAT header specification (Affymetrix Data File Formats, April
  # 2006), we know that the date and the timestamp is 18 characters long.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   nTemp <- 7;
##   nPower <- 4;
##   nTimestamp <- 18;
##   # Expected start position
##   pos <- nTemp + 1 + nPower + 1;
##   # ...however, the files we have start at the next position. /HB 2006-12-01
##   pos <- pos + 1;
##   timestamp <- substring(header, first=pos, last=pos+nTimestamp-1);

  timestamp <- trim(timestamp); # Unnecessary?
  res <- strptime(timestamp, format=format, ...);
  attr(res, "text") <- timestamp;
  res;
}, private=TRUE)



setMethodS3("nbrOfCells", "AffymetrixCelFile", function(this, ...) {
  getHeader(this)$total;
})


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
}, private=TRUE)



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
# \section{Caching}{
#   CEL data is neither cached in memory nor on file by this method.
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
setMethodS3("readUnits", "AffymetrixCelFile", function(this, units=NULL, cdf=NULL, ..., stratifyBy=NULL, force=FALSE, cache=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve CDF structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cdf)) {
    suppressWarnings({
      cdf <- readUnits(getCdf(this), units=units, stratifyBy=stratifyBy);
    });
  } 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    res <- readCelUnits(this$.pathname, cdf=cdf, ...);
  })

  res;
}, private=TRUE)


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
#   \item{...}{Not used.}
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
}, private=TRUE)



###########################################################################/**
# @RdocMethod clearData
#
# @title "Clears all or a subset of the fields in a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fields}{A @character @vector of fields to be cleared.}
#   \item{value}{A @numeric value to be written over the data.}
#   \item{...}{Additional arguments passed to the
#      @see "affxparser::updateCelUnits" methods.}
#   \item{.forSure}{If not @TRUE, an exception is thrown asking if the
#      method was called by mistake.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) the names of the fields cleared.
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
setMethodS3("clearData", "AffymetrixCelFile", function(this, fields=c("intensities", "stdvs", "pixels"), value=0, ..., .forSure=FALSE, verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  if (!identical(.forSure, TRUE))
    throw("Did you call clearData() by mistake? If not, use .forSure=TRUE.");

  # Nothing do to?
  if (length(fields) == 0) {
    verbose && cat(verbose, "No fields to be cleared.");
    return(invisible(fields));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clear
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Clearing Affymetrix CEL file");
  verbose && cat(verbose, "Fields to be cleared: ", paste(fields, collapse=", "));
  bfr <- rep(value, length.out=nbrOfCells(this));
  intensities <- stdvs <- pixels <- NULL;
  if ("intensities" %in% fields)
    intensities <- bfr;
  if ("stdvs" %in% fields)
    stdvs <- bfr;
  if ("pixels" %in% fields)
    pixels <- bfr;
  updateCel(getPathname(this), intensities=bfr, stdvs=bfr, pixels=bfr);
  verbose && exit(verbose);

  invisible(fields);
}, static=TRUE, private=TRUE)




setMethodS3("[", "AffymetrixCelFile", function(this, units=NULL, drop=FALSE) {
  data <- readUnits(this, units=units);
  if (drop && length(data) == 1)
    data <- data[[1]];
  data;
})

setMethodS3("[[", "AffymetrixCelFile", function(this, unit=NULL) {
  this[units=unit, drop=TRUE];
})



###########################################################################/**
# @RdocMethod getData
#
# @title "Gets all or a subset of the fields in a CEL file"
#
# \description{
#  @get "title" for all or a subset of the cells.
# }
#
# @synopsis
#
# \arguments{
#   \item{indices}{A @numeric @vector of cell indices.  If @NULL, all cells
#     are considered.}
#   \item{fields}{Names of fields to be retrieved.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCel".}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @data.frame of the fields requested.
# }
#
# \section{Caching}{
#   Neither in-memory nor on-file caching is done by this method.
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
setMethodS3("getData", "AffymetrixCelFile", function(this, indices=NULL, fields=c("xy", "intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
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

  # Argument 'fields':
  fields <- intersect(fields, c("xy", "intensities", "stdvs", "pixels"));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cVerbose <- -(as.numeric(verbose) + 50);
  args <- list(
    filename=this$.pathname, indices=indices, 
    readHeader=FALSE, 
    readIntensities=("intensities" %in% fields), 
    readStdvs=("stdvs" %in% fields), 
    readPixels=("pixels" %in% fields), 
    readXY=("xy" %in% fields), 
    readOutliers=FALSE, 
    readMasked=FALSE,
    ...,
    verbose=cVerbose
  );
  fcn <- get("readCel", mode="function");
  keep <- intersect(names(args), names(formals(fcn)));
  args <- args[keep];
  cel <- do.call("readCel", args=args);
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split (x,y)
  isXY <- which(fields == "xy");
  if (length(isXY) > 0) {
    fields <- as.list(fields);
    fields[[isXY]] <- c("x", "y");
    fields <- unlist(fields, use.names=TRUE);
  }

  # Keep only requested fields
  if (!identical(names(cel), fields)) {
    cel <- cel[fields];
  }

  # Return as data frame
  attr(cel, "row.names") <- seq_len(length(cel[[1]]));
  class(cel) <- "data.frame";

  cel;
}, private=TRUE)


setMethodS3("range", "AffymetrixCelFile", function(this, ..., na.rm=TRUE) {
  x <- getData(this, ...);
  range(x, na.rm=na.rm);
}, private=TRUE)


setMethodS3("getRectangle", "AffymetrixCelFile", function(this, xrange=c(0,Inf), yrange=c(0,Inf), fields=c("intensities", "stdvs", "pixels"), ...) {
  readCelRectangle(this$.pathname, xrange=xrange, yrange=yrange, readIntensities=("intensities" %in% fields), readStdvs=("stdvs" %in% fields), readPixels=("pixels" %in% fields));
}, private=TRUE)

setMethodS3("setAttributeXY", "AffymetrixCelFile", function(this, value, ...) {
  # Argument 'value':
  if (is.null(value)) {
    # Nothing todo?
    return();
  }

  pattern <- "^(X*)(Y*)$";
  if (regexpr(pattern, value) == -1) {
    throw("The value of argument 'value' is unrecognized: ", value);
  }

  # Parse and count
  n23 <- gsub(pattern, "\\1", value);
  n24 <- gsub(pattern, "\\2", value);
  n23 <- nchar(n23);
  n24 <- nchar(n24);
  setAttributes(this, n23=n23, n24=n24);
})

setMethodS3("getAttributeXY", "AffymetrixCelFile", function(this, ...) {
  n23 <- getAttribute(this, "n23", 0);
  n24 <- getAttribute(this, "n24", 0);
  xyTag <- paste(c(rep("X", n23), rep("Y", n24)), collapse="");
  xyTag;
})

setMethodS3("hasAttributeXY", "AffymetrixCelFile", function(this, values, ...) {
   xyTag <- getAttributeXY(this);
   (xyTag %in% values);
})


setMethodS3("parseTagsAsAttributes", "AffymetrixCelFile", function(this, ...) {
  newAttrs <- NextMethod("parseTagsAsAttributes", this, ...);

  # Get all XY, XX, XXX etc tags
  tags <- getTags(this, pattern="^X*Y*$");
  newAttrs <- c(newAttrs, setAttributeXY(this, tags));

  # Return nothing
  invisible(newAttrs);
})



############################################################################
# HISTORY:
# 2007-03-05
# o Added parseTagsAsAttributes().
# o Added setAttributeXY(), getAttributeXY(), and hasAttributeXY().
# 2007-02-12
# o Now getData() is using do.call() because it is faster. Unused arguments
#   are still ignored.
# 2007-02-04
# o Now getData() is call readCel() using doCall() so that unused arguments
#   in '...' are ignored.
# 2007-02-03
# o BUG FIX: getTimestamp() assumed a fix location in the CEL v3 header,
#   but that did not work for dChip exported CEL files.  Now, a US date
#   pattern is assumed and searched for.
# 2007-01-12 /KS
# o Moved image270() and writeSpatial() to AffymetrixCelFile.PLOT.R.
# 2006-12-18 /KS
# o Add "takeLog" argument (logical) to image270.  If true, take the log2
#   before plotting.  Can be more informative than natural scale.
# 2006-12-14
# o Removed getSampleName() which gives the same as getName().
# 2006-12-11
# o Now the timestamp is also reported for singel CEL files.
# o BUG FIX: getHeaderV3() would throw an error if there was an empty V3 
#   header fields.  This was the reason why getTimestamp() gave an error
#   on some 100K chips.
# 2006-12-01
# o Added getTimestamp().
# 2006-11-28
# o Arguments 'force' and 'cache' has to be in readUnits() to avoid being
#   passed from calls of subclasses.
# 2006-10-23
# o Update default value for argument 'fields' in getData().
# 2006-10-22
# o In order to speed up fromFile(), the CEL header is not read anymore.
# 2006-10-06
# o make sure cdf association is inherited
# 2006-08-28
# o Renamed getFields() to getData() because getFields() is "reserved"
#   for use in the Object class.
# 2006-08-27
# o Added nbrOfCells() because it is so common.
# o Added createFrom() which utilizes new functions copyFile() and 
#   clearData(). It is also no longer static. This is more generic and 
#   cleaner.  The new clearData() does also not require the CDF file 
#   (in case that should be missing).
# 2006-08-25
# o Renamed getIntensities() to getFields() which returns a data frame.
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
