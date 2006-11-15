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
#   \item{encodeFunction}{A @function taking a single @list structure
#      as its argument.}
#   \item{decodeFunction}{A @function taking a single @list structure
#      as its argument.}
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
setConstructorS3("ParameterCelFile", function(..., encodeFunction=NULL, decodeFunction=NULL) {
  extend(AffymetrixCelFile(...), "ParameterCelFile",
    "cached:.readUnitsCache" = NULL,
    encodeFunction = encodeFunction,
    decodeFunction = decodeFunction
  )
})


setMethodS3("setEncodeFunction", "ParameterCelFile", function(this, fcn, ...) {
  if (is.null(fcn)) {
  } else if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", mode(fcn));
  }
  this$encodeFunction <- fcn;
  invisible(this);
})

setMethodS3("setDecodeFunction", "ParameterCelFile", function(this, fcn, ...) {
  if (is.null(fcn)) {
  } else if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", mode(fcn));
  }
  this$decodeFunction <- fcn;
  invisible(this);
})

# There was a lot of overhead for calling functions in the previous 
# encoding/decoding mechanism where encode() called encodeUnit() for
# every unit individually. Same for decode() and decodeUnits(). The
# new mechanism skips the uncodeUnit() and decodeUnit() step.
setMethodS3("encode", "ParameterCelFile", function(this, units, ...) {
  encodeUnitGroup <- this$encodeFunction;
  if (!is.null(encodeUnitGroup)) {
    units <- lapply(units, FUN=lapply, encodeUnitGroup, ...);
  }
  units;
}, protected=TRUE)

setMethodS3("decode", "ParameterCelFile", function(this, units, ...) {
  decodeUnitGroup <- this$decodeFunction;
  if (!is.null(decodeUnitGroup)) {
    units <- lapply(units, FUN=lapply, decodeUnitGroup, ...);
  }
  units;
}, protected=TRUE)


# Previous solutions:
# setMethodS3("encodeUnit", "ParameterCelFile", function(this, unit, ...) {
#   encodeUnitGroup <- this$encodeFunction;
#   if (!is.null(encodeUnitGroup)) {
#     unit <- lapply(unit, FUN=encodeUnitGroup, ...);
#   }
#   unit;
# }, protected=TRUE)
# 
# setMethodS3("encode", "ParameterCelFile", function(this, units, ...) {
#   if (!is.null(this$encodeFunction)) {
#     units <- lapply(units, FUN=function(unit) encodeUnit(this, unit, ...));
#   }
#   units;
# }, protected=TRUE)
#
# setMethodS3("decodeUnit", "ParameterCelFile", function(this, unit, ..., verbose=FALSE) {
#   decodeUnitGroup <- this$decodeFunction;
#   if (!is.null(decodeUnitGroup)) {
#     unit <- lapply(unit, FUN=decodeUnitGroup, ...);
#   }
#   unit;
# }, protected=TRUE)
# 
# setMethodS3("decode", "ParameterCelFile", function(this, units, ..., verbose=FALSE) {
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Validate arguments
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Argument 'verbose':
#   verbose <- Arguments$getVerbose(verbose);
# 
# 
#   verbose && enter(verbose, "Decoding ", length(units), " units");
# 
#   if (!is.null(this$decodeFunction)) {
#     # Clever decoding... If three units across the file are not modified
#     # by the decoding method, don't decode the others either.
#     someUnits <- units[seq(from=1,to=length(units),length=3)];
#     someUnits2 <- lapply(someUnits, FUN=function(unit) decodeUnit(this, unit, ...));
#     if (!identical(someUnits, someUnits2)) {
#       units <- lapply(units, FUN=function(unit) decodeUnit(this, unit, ...));
#     } else {
#       verbose && cat(verbose, "Decoding not needed. No changes detected.");
#     }
#   }
#   verbose && exit(verbose);
# 
#   units;
# }, protected=TRUE)



setMethodS3("readUnits", "ParameterCelFile", function(this, ..., readStdvs=FALSE, readPixels=FALSE, stratifyBy=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating hashcode key for cache");
  key <- digest(list(..., readStdvs=readStdvs, readPixels=readPixels, stratifyBy=stratifyBy));
  verbose && exit(verbose);
  if (!force) {
    verbose && enter(verbose, "Trying to obtain cached data");
    res <- this$.readUnitsCache[[key]];
    verbose && exit(verbose);
    if (!is.null(res)) {
      verbose && cat(verbose, "readUnits.ParameterCelFile(): Returning cached data");
      return(res);
    }
  }
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve and decoding data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading units");
  units <- NextMethod("readUnits", this, ..., readStdvs=readStdvs, readPixels=readPixels, stratifyBy=stratifyBy, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Decoding ", length(units), " units");
  units <- decode(this, units, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "readUnits.ParameterCelFile(): Updating cache");
  this$.readUnitsCache <- list();
  this$.readUnitsCache[[key]] <- units;

  units;
});


setMethodS3("updateUnits", "ParameterCelFile", function(this, data, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating units");
 
  verbose && enter(verbose, "Encoding units");
  data <- encode(this, data);
  verbose && exit(verbose);

  NextMethod("updateUnits", this, cdf=cdf, data=data, ..., verbose=less(verbose));

  verbose && exit(verbose);

  invisible(data);
}, protected=TRUE);





############################################################################
# HISTORY:
# 2006-11-14
# o Removed encode- and decodeUnit(). It caused a substantial overhead;
#   it is about 2-3 faster letting encode() call the encode function
#   directly. I made sure that encode() and decode() gives the same results
#   before and after the update. /HB
# o Added 'verbose' to updateUnits().
# 2006-09-10
# o Now the encode and decode functions for a unit group are made into
#   fields of this class.  This way we don't have to create a special
#   ParameterCelFile class for each kind of model.
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
