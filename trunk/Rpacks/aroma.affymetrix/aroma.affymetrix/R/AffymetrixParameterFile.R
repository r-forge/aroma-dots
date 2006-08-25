###########################################################################/**
# @RdocClass AffymetrixParameterFile
#
# @title "The AffymetrixParameterFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixParameterFile object represents parameter estimates for
#  a single Affymetrix array.
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
#*/###########################################################################
setConstructorS3("AffymetrixParameterFile", function(...) {
  extend(AffymetrixFile(...), "AffymetrixParameterFile",
    .cdf = NULL,
    .fields = list()
  )
})


###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixParameterFile object from a file"
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
#  Returns an @see "AffymetrixParameterFile" object.  
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
setMethodS3("fromFile", "AffymetrixParameterFile", function(static, filename, path=NULL, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  isApf <- (regexpr("[.](a|a)(p|P)(f|F)$", pathname) != -1);
  if (!isCel) {
    throw("Could not read APF file. File format error: ", pathname);
  }

  # Try to read the header assuming
  header <- readApfHeader(pathname);

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("readHeader", "AffymetrixParameterFile", function(this, ...) {
  readInteger <- function(con, n=1, ...) {
    readBin(con, what=integer(), n=n, size=4, signed=TRUE, endian="little");
  } # readInteger()

  # Open file
  con <- file(this$.pathname, open="rb");
  on.exit(close(con));

  header <- list();

  nchars <- readInteger(con);
  header$comment <- readChar(con, nchars=nchars);

  nbrOfFields <- readInteger(con);
  names <- vector("character", length=nbrOfFields);
  types <- vector("character", length=nbrOfFields);
  whats <- vector("character", length=nbrOfFields);
  sizes <- vector("integer", length=nbrOfFields);
  lengths <- vector("integer", length=nbrOfFields);

  for (kk in seq(length=nbrOfFields)) {
    nchars <- readInteger(con);
    names[kk] <- readChar(con, nchars=nchars);
    nchars <- readInteger(con);
    types[kk] <- readChar(con, nchars=nchars);
    res <- whatDataType(types[kk]);
    whats[kk] <- res$what;
    sizes[kk] <- res$size;
    lengths[kk] <- readInteger(con);
  }

  dataOffset <- seek(con, rw="read");
  offsets <- dataOffset + (cusum(sizes*lengths) - sizes[1]*lengths[1]);

  names(types) <- names;
  names(whats) <- names;
  names(sizes) <- names;
  names(lengths) <- names;
  names(offsets) <- offsets;

  header$names <- names;
  header$types <- types;
  header$whats <- whats;
  header$sizes <- sizes;
  header$lengths <- lengths;
  header$offsets <- offsets;

  header;
}, protected=TRUE)


setMethodS3("readField", "AffymetrixParameterFile", function(this, name, subset=NULL, ...) {
  header <- getHeader(this);

  # Argument 'name':
  if (name %in% header$names)
    throw("No such field: ", name);

  # Argument 'subset':
  length <- header$lengths[name];
  if (!is.null(subset)) {
    if (any(subset < 1 | subset > length))
      throw(sprintf("Argument 'subset' is out of range [1,%d].", length));
  }

  # Open file
  con <- file(this$.pathname, open="rb");
  on.exit(close(con));

  # Number of elements to read
  if (is.null(subset)) {
    n <- length;
  } else {
    n <- length(subset);
  }

  # Data type to read
  what <- header$whats[name];

  # Nothing more to do?
  if (n == 0)
    return(vector(what, 0));

  # Move the reading file position
  offset <- header$offsets[name];
  seek(con, where=offset, origin="start", rw="read");

  # Size of data type to be read
  size <- header$sizes[name];

  # Read everything?
  if (is.null(subset)) {
    x <- readBin(con, what=what, size=size, endian="little", n=n);
    return(x);
  }

  # Read a subset:
  # First, read it in an optimal order
  o <- order(subset);

  # Ordered subset
  subset <- subset(o);

  # Skip to the first element to be read
  seek(con, where=subset[1], origin="current", rw="read");
  subset <- subset - subset[1];

  # Read one contiguous block
  x <- readBin(con, what=what, size=size, endian="little", n=max(subset));

  # Extract subset of interest
  x <- x[subset + 1];

  # Re-reorder data
  x <- x[order(o)];

  x;
}, protected=TRUE)

setMethodS3("getHeader", "AffymetrixParameterFile", function(this, ...) {
  header <- this$.header;
  if (is.null(header))
    header <- this$.header <- readHeader(this);
  header;
}, protected=TRUE)

setMethodS3("getChipType", "AffymetrixParameterFile", function(this, ...) {
  getHeader(this)$chiptype;
}, protected=TRUE)

setMethodS3("getCdf", "AffymetrixParameterFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    cdf <- this$.cdf <- AffymetrixCdfFile$fromChipType(getChipType(this));
  }
  cdf;
})

setMethodS3("getIntensities", "AffymetrixParameterFile", function(this, indices=NULL, ..., verbose=FALSE) {
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


setMethodS3("getUnits", "AffymetrixParameterFile", function(this, ...) {
  readCelUnits(this$.pathname, cdf=getPathname(getCdf(this)), ...);
}, protected=TRUE);


############################################################################
# HISTORY:
# 2006-08-19
# o Created.
############################################################################
