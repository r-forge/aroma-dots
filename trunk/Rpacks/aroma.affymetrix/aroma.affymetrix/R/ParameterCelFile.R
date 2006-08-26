###########################################################################/**
# @RdocClass ParameterCelFile
#
# @title "The ParameterCelFile class"
#
# \description{
#  @classhierarchy
#
#  A ParameterCelFile object represents parameter estimates.
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
# \section{File format}{
#   The idea behind this class is store data fields which by nature have 
#   one value per probe (per field) in CEL files.  A perfect example is to
#   store probe-affinity estimates and their standard deviations.  There
#   is one probe affinity per probe so the structure of a CEL file (and
#   its coupled CDF file) is well suited to read/write such information.
#   
#   Consider a unit group with L probes.  A CEL file stores 
#   \code{intensities} (L floats), \code{stdvs} (L floats), and 
#   \code{pixels} (L integers).  Thus, for each probe l=1,...,L, a
#   (float, float, integer) tuple is stored.  We can use this for any
#   information we want.  If we want a slightly different structure,
#   we can choose to encode/decode our structure/information to fit the
#   structure of the CEL file.  This abstract class provides transparent 
#   methods for encoding and decoding such information through methods
#   @seemethod "encodeUnitGroup" and @seemethod "decodeUnitGroup".
#   By subclassing you can implement different types of data structures.
# }
#
# @author
#
# @keyword "IO"
#*/###########################################################################
setConstructorS3("ParameterCelFile", function(...) {
  extend(AffymetrixCelFile(...), "ParameterCelFile")
}, abstract=TRUE)




setMethodS3("createFrom", "ParameterCelFile", function(static, dataSet, filename, path=NULL, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", 
                                                      class(dataSet)[1]);
    }
    if (length(dataSet) == 0) {
      throw("Cannot create ", class(static)[1], " object. The ", class(dataSet), " does not contain any files.");
    }
  }

  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Don't create if already there
  if (isFile(pathname))
    return(pathname);

  verbose && enter(verbose, "Creating file for probe-affinity estimates.");
  # 1. Get the pathname to get first CEL file
  copyFrom <- getPathname(dataSet[[1]]);
  file.copy(copyFrom, pathname);

  # 2. Clear the file
  cdf <- getCdf(dataSet);
  bfr <- double(nbrOfCells(cdf));
  updateCel(pathname, intensities=bfr, stdvs=bfr, pixels=bfr);

  verbose && exit(verbose);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("encodeUnitGroup", "ParameterCelFile", abstract=TRUE, static=TRUE, protected=TRUE);

setMethodS3("encodeUnit", "ParameterCelFile", function(static, unit, ...) {
  lapply(unit, FUN=function(group) encodeUnitGroup(static, group, ...));
}, protected=TRUE)

setMethodS3("encode", "ParameterCelFile", function(static, units, ...) {
  lapply(units, FUN=function(unit) encodeUnit(static, unit, ...));
}, protected=TRUE)

setMethodS3("decodeUnitGroup", "ParameterCelFile", function(static, intensities=NULL, stdvs=NULL, pixels=NULL, ...) {
  list(intensities=intensities, stdvs=stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)

setMethodS3("decodeUnit", "ParameterCelFile", function(static, unit, ...) {
  lapply(unit, FUN=function(group) decodeUnitGroup(static, group, ...));
}, protected=TRUE)

setMethodS3("decode", "ParameterCelFile", function(static, units, ...) {
  lapply(units, FUN=function(unit) decodeUnit(static, unit, ...));
}, protected=TRUE)


setMethodS3("readUnits", "ParameterCelFile", function(this, ..., readStdvs=TRUE, readPixels=TRUE, stratifyBy=NULL) {
  units <- NextMethod("readUnits", this, ..., readStdvs=readStdvs, readPixels=readPixels, stratifyBy=stratifyBy);
  decode(this, units);
});

setMethodS3("updateUnits", "ParameterCelFile", function(this, data, ...) {
  data <- encode(this, data);
  NextMethod("updateUnits", this, data=data, ...);
  invisible(data);
}, protected=TRUE);





############################################################################
# HISTORY:
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added updateUnits().
# 2006-08-21
# o Created.
############################################################################
