###########################################################################/**
# @RdocClass AromaSignalBinaryFile
#
# @title "The AromaSignalBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaSignalBinaryFile is a @see "AromaTabularBinaryFile".
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
#   @see "AromaTabularBinaryFile".
# }
#*/########################################################################### 
setConstructorS3("AromaSignalBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), c("AromaSignalBinaryFile",
                                              uses("AromaPlatformInterface")),
    "cached:.unf" = NULL,
    "cached:.ugp" = NULL
  );
})




setMethodS3("as.character", "AromaSignalBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));

  class(s) <- class;
  s;
})


setMethodS3("fromFile", "AromaSignalBinaryFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
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



setMethodS3("getFilenameExtension", "AromaSignalBinaryFile", function(static, ...) {
  "asb";
}, static=TRUE, protected=TRUE)



setMethodS3("nbrOfUnits", "AromaSignalBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("allocate", "AromaSignalBinaryFile", function(static, ..., platform, chipType, types="double", sizes=4, signed=TRUE, footer=list()) {
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Create tabular binary file
  res <- allocate.AromaTabularBinaryFile(static, generic="allocate", ...,
                                  types=types, sizes=sizes, signed=signed);


  # Write attributes to footer
  attrs <- list(
    platform=platform, 
    chipType=chipType
  );
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getPlatform", "AromaSignalBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  res <- footer$platform;

  if (!is.null(res)) {
    res <- as.character(res);
    res <- unlist(strsplit(res, split="[\t]"));
    res <- trim(res);
  }

  res;
})

setMethodS3("getChipType", "AromaSignalBinaryFile", function(this, ...) {
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


setMethodS3("getUnitNamesFile", "AromaSignalBinaryFile", function(this, force=FALSE, ...) {
  unf <- this$.unf;
  if (force || is.null(unf)) {
    platform <- getAromaPlatform(this);
    unf <- getUnitNamesFile(platform);
    this$.unf <- unf;
  }

  unf;
})


setMethodS3("getAromaUgpFile", "AromaSignalBinaryFile", function(this, ..., force=FALSE) {
  ugp <- this$.ugp;
  if (force || is.null(ugp)) {
    chipType <- getChipType(this, ...);
    ugp <- AromaUgpFile$byChipType(chipType);
    this$.ugp <- ugp;
  }
  ugp;
}) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
