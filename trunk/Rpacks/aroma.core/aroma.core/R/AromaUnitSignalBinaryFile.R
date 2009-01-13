###########################################################################/**
# @RdocClass AromaUnitSignalBinaryFile
#
# @title "The AromaUnitSignalBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitSignalBinaryFile is a @see "AromaTabularBinaryFile".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   @see "aroma.core::AromaTabularBinaryFile".
# }
#*/########################################################################### 
setConstructorS3("AromaUnitSignalBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), c("AromaUnitSignalBinaryFile",
                                              uses("AromaPlatformInterface")),
    "cached:.unf" = NULL,
    "cached:.ugp" = NULL
  );
})




setMethodS3("as.character", "AromaUnitSignalBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));

  class(s) <- class;
  s;
})


setMethodS3("fromFile", "AromaUnitSignalBinaryFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
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

  res <- newInstance(static, filename=pathname, path=NULL, ...);

  res;
})



setMethodS3("getFilenameExtension", "AromaUnitSignalBinaryFile", function(static, ...) {
  "asb";
}, static=TRUE, protected=TRUE)



setMethodS3("nbrOfUnits", "AromaUnitSignalBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("allocate", "AromaUnitSignalBinaryFile", function(static, ..., platform, chipType, types="double", sizes=4, signed=TRUE, footer=list()) {
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Create tabular binary file
  res <- allocate.AromaTabularBinaryFile(static, generic="allocate", ...,
                                  types=types, sizes=sizes, signed=signed);


  # Write attributes to footer
  attrs <- list(
    createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
    platform=platform, 
    chipType=chipType
  );
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE)



setMethodS3("extractMatrix", "AromaUnitSignalBinaryFile", function(this, units=NULL, rows=units, ...) {
  NextMethod("extractMatrix", rows=units, ...);  
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getPlatform", "AromaUnitSignalBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  res <- footer$platform;

  if (!is.null(res)) {
    res <- as.character(res);
    res <- unlist(strsplit(res, split="[\t]"));
    res <- trim(res);
  }

  res;
})

setMethodS3("getChipType", "AromaUnitSignalBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  res <- footer$chipType;

  if (is.null(res)) {
    throw("File format error: This ", class(this)[1], " file does not contain information on chip type in the file footer: ", getPathname(this));
  }

  res <- as.character(res);
  res <- unlist(strsplit(res, split="[\t]"));
  res <- trim(res);

  res;
})


setMethodS3("getUnitNamesFile", "AromaUnitSignalBinaryFile", function(this, force=FALSE, ...) {
  unf <- this$.unf;
  if (force || is.null(unf)) {
    platform <- getAromaPlatform(this);
    unf <- getUnitNamesFile(platform);
    this$.unf <- unf;
  }

  unf;
})



setMethodS3("allocateFromUnitNamesFile", "AromaUnitSignalBinaryFile", function(static, unf, ...) {
  # Argument 'unf':
  className <- "UnitNamesFile";
  if (!inherits(unf, className)) {
    throw("Argument 'unf' is not of class ", className, ": ", class(unf)[1]);
  }

  platform <- getPlatform(unf);
  chipType <- getChipType(unf);
  nbrOfRows <- nbrOfUnits(unf);
  
  allocate.AromaUnitSignalBinaryFile(static, ..., nbrOfRows=nbrOfRows, platform=platform, chipType=chipType);
}, static=TRUE)


setMethodS3("getAromaUgpFile", "AromaUnitSignalBinaryFile", function(this, ..., validate=FALSE, force=FALSE) {
  ugp <- this$.ugp;
  if (force || is.null(ugp)) {
    chipType <- getChipType(this, ...);
    ugp <- AromaUgpFile$byChipType(chipType, validate=validate);
    this$.ugp <- ugp;
  }
  ugp;
}) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2009-01-12
# o Added extractMatrix() accepting argument 'units'.  This will then
#   also work for the corresponding set of files.
# 2009-01-05
# o Renamed from AromaSignalBinaryFile to AromaUnitSignalBinaryFile.
# 2008-05-24
# o Now allocate() of AromaSignalBinaryFile adds footer 'createdOn'.
# 2008-05-11
# o Created.
############################################################################
