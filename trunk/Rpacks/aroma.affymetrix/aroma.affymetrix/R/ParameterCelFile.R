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




###########################################################################/**
# @RdocMethod encodeUnitGroup
#
# @title "Returns a unit group with fields as in a CEL file"
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
#  Returns a named @list structure with elements \code{intensities},
#  \code{stdvs}, and \code{pixels}, all of equal length.
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
setMethodS3("encodeUnitGroup", "ParameterCelFile", abstract=TRUE, static=TRUE, protected=TRUE);




setMethodS3("encodeUnit", "ParameterCelFile", function(static, unit, ...) {
  lapply(unit, FUN=function(group) encodeUnitGroup(static, group, ...));
}, protected=TRUE)

setMethodS3("encode", "ParameterCelFile", function(static, units, ...) {
  lapply(units, FUN=function(unit) encodeUnit(static, unit, ...));
}, protected=TRUE)

setMethodS3("decodeUnitGroup", "ParameterCelFile", function(static, intensities=NULL, stdvs=NULL, pixels=NULL, ...) {
  res <- list();
  if (!is.null(intensities))
    res$intensities <- intensities;
  if (!is.null(stdvs))
    res$stdvs <- stdvs;
  if (!is.null(pixels))
    res$pixels <- pixels;
  res;
}, static=TRUE, protected=TRUE)

setMethodS3("decodeUnit", "ParameterCelFile", function(static, unit, ...) {
  lapply(unit, FUN=function(group) decodeUnitGroup(static, group, ...));
}, protected=TRUE)

setMethodS3("decode", "ParameterCelFile", function(static, units, ...) {
  lapply(units, FUN=function(unit) decodeUnit(static, unit, ...));
}, protected=TRUE)


setMethodS3("readUnits", "ParameterCelFile", function(this, ..., readStdvs=FALSE, readPixels=FALSE, stratifyBy=NULL) {
  units <- NextMethod("readUnits", this, ..., readStdvs=readStdvs, readPixels=readPixels, stratifyBy=stratifyBy);
  decode(this, units);
});

setMethodS3("updateUnits", "ParameterCelFile", function(this, data, cdf=NULL, ...) {
  data <- encode(this, data);
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
  invisible(data);
}, protected=TRUE);





############################################################################
# HISTORY:
# 2006-08-27
# o Moved createFrom() to the AffymetrixCelFile class.
# 2006-08-26
# o Now createFrom() takes a CEL file and not a CEL set.
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added updateUnits().
# 2006-08-21
# o Created.
############################################################################
