###########################################################################/**
# @RdocClass AromaUnitTabularBinaryFile
#
# @title "The AromaUnitTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  A AromaUnitTabularBinaryFile is an @see "AromaTabularBinaryFile" with
#  the constraint that the rows map one-to-one to, and in the same order as,
#  the units in a annotation chip type file (e.g. CDF file).  
#  The (full) chip type of the annotation chip type file is given by the
#  mandatory file footer \code{chipType}.
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
setConstructorS3("AromaUnitTabularBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), c("AromaUnitTabularBinaryFile",
                                              uses("AromaPlatformInterface")),
    "cached:.unf" = NULL
  );
})


setMethodS3("as.character", "AromaUnitTabularBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  n <- length(s);
  s <- s[c(1:(n-2), n, n-1)];

  class(s) <- class;
  s;
})



setMethodS3("clearCache", "AromaUnitTabularBinaryFile", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".unf")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...); 
}, private=TRUE)


setMethodS3("getFilenameExtension", "AromaUnitTabularBinaryFile", static=TRUE, abstract=TRUE);

setMethodS3("nbrOfUnits", "AromaUnitTabularBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("getPlatform", "AromaUnitTabularBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  platform <- footer$platform;

  if (is.null(platform)) {
    # AD HOC: If there is no platform information in the file, then assume
    # it is an Affymetrix platform.  Newer files are allocated with 
    # platform information.  Older files have only be created for the 
    # Affymetrix platform. /HB 2008-05-18.
    platform <- "Affymetrix";
    warning(sprintf("%s does not have a 'platform' footer attribute. Assuming platform 'Affymetrix': %s", class(this)[1], getPathname(this)));
  } 

  if (!is.null(platform)) {
    platform <- as.character(platform);
    platform <- unlist(strsplit(platform, split="[\t]"));
    platform <- trim(platform);
  }

  platform;
})


setMethodS3("getChipType", "AromaUnitTabularBinaryFile", function(this, fullname=TRUE, .old=FALSE, ...) {
  footer <- readFooter(this);
  chipType <- footer$chipType;

  if (is.null(chipType) && !.old) {
    throw("File format error: This ", class(this)[1], " file does not contain information on chip type in the file footer.  This is because the file is of an older file format an is no longer supported.  Please update to a more recent version: ", getPathname(this));
  }

  if (.old) {
    # Keep backward compatible for a while. /HB 2008-01-19
    if (fullname) {
      # Currently these files do not contain fullname chiptype information.
      # Instead we have to search for a match CDF and use its chip type.
      # Note, the CDF returned will depend of which CDF exists in the search
      # path, if at all.  AD HOC! /HB 2007-12-10
      unf <- getUnitNamesFile(this, ...);
      chipType <- getChipType(unf, fullname=fullname);
    } else {
      chipType <- getName(this, ...);
    }
    chipType <- trim(chipType);
    if (nchar(chipType) == 0)
      throw("File format error: The inferred chip type is empty.");
  } else {
    if (!fullname) {
      chipType <- gsub(",.*", "", chipType);
    }
    chipType <- trim(chipType);
    if (nchar(chipType) == 0)
      throw("File format error: The chip type according to the file footer is empty.");
  }

  if (!is.null(chipType)) {
    chipType <- as.character(chipType);
    chipType <- unlist(strsplit(chipType, split="[\t]"));
    chipType <- trim(chipType);
  }

  chipType;
})



setMethodS3("getUnitNamesFile", "AromaUnitTabularBinaryFile", function(this, force=FALSE, ...) {
  unf <- this$.unf;
  if (force || is.null(unf)) {
    platform <- getAromaPlatform(this, ...);
    unf <- getUnitNamesFile(platform);
    this$.unf <- unf;
  }

  unf;
})


setMethodS3("byChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, validate=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Locating ", class(static)[1]);

  pathname <- findByChipType(static, chipType=chipType, tags=tags, 
                                                     firstOnly=TRUE, ...);
  if (is.null(pathname)) {
    throw("Could not locate a file for this chip type: ", 
                                   paste(c(chipType, tags), collapse=","));
  }

  verbose && cat(verbose, "Located file: ", pathname);

  # Create object
  res <- newInstance(static, pathname);

  verbose && exit(verbose);

  res;
}, static=TRUE)



setMethodS3("fromChipType", "AromaUnitTabularBinaryFile", function(static, ...) {
  byChipType(static, ...);
}, deprecated=TRUE)


setMethodS3("findByChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  # Strip 'monocell' parts
  parts <- parts[parts != "monocell"];
  chipType <- parts[1];
  tags <- parts[-1];
  fullname <- paste(c(chipType, tags), collapse=",");

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
#  args$firstOnly <- FALSE;
#  str(args);
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s.*[.]%s[.]lnk$", chipType, ext);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)




setMethodS3("indexOfUnits", "AromaUnitTabularBinaryFile", function(this, names, ...) {
  # Map the unit names to the ones in the unit names file
  unf <- getUnitNamesFile(this);
  unitNames <- getUnitNames(unf);
  idxs <- match(names, unitNames);
  idxs;
}, protected=TRUE)




setMethodS3("allocateFromUnitNamesFile", "AromaUnitTabularBinaryFile", function(static, unf, path=getPath(unf), tags=NULL, footer=list(), ...) {
  # Argument 'unf':
  if (!inherits(unf, "UnitNamesFile")) {
    throw("Argument 'unf' is not an UnitNamesFile: ", class(unf)[1]);
  }

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(unf);

  # Exclude 'monocell' tags (AD HOC)
  chipType <- gsub(",monocell", "", chipType);

  # Get platform
  platform <- getPlatform(unf);

  # Number of units
  unf <- nbrOfUnits(unf);

  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create tabular binary file
  res <- allocate(static, filename=filename, path=path, 
                                                nbrOfRows=nbrOfUnits, ...);


  # Write attributes to footer
  attrs <- list(platform=platform, chipType=chipType);
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE)



############################################################################
# HISTORY:
# 2008-05-19
# o Added getPlatform().
# o Added platform-independent allocateFromUnitNamesFile() which now also
#   writes footer attribute 'platform'.
# 2008-02-13
# o Added and updated Rdoc comments.
# 2008-01-19
# o Now AromaUnitTabularBinaryFile gets the chip type from the file footer.
# o ROBUSTNESS: Now fromChipType() of AromaUnitTabularBinaryFile validates
#   that the number of units in the located file match the number of units
#   in the CDF located using the same search parameters.
# 2007-12-10
# o Currently a AromaUnitTabularBinaryFile (e.g. AromaUgpFile) does not
#   contain information about the "fullname" chip type, but only the basic
#   chip-type name, e.g. we cannot infer the full chip-type name from 
#   'GenomeWideSNP_5,Full,r2.ugp', but only 'GenomeWideSNP_5'. The fullname
#   should be the same as the full chip-type name of the CDF used to define
#   the the unit map, e.g. 'GenomeWideSNP_5,Full.CDF'.
#   We should add a header (or footer) field in the file format that 
#   indicates the full chip type.  
#   However, until that is done, the best we can do is to turn to the ad
#   hoc solution of scanning for the CDF with the longest matching fullname,
#   if both 'GenomeWideSNP_5,Full.CDF' and 'GenomeWideSNP_5.CDF' exists,
#   the we match the former to 'GenomeWideSNP_5,Full,r2.ugp'.  The fullname
#   chip type of the UGP is then full chip-type name of the CDF.  NOTE,
#   there is major drawback with this.  If the user deletes the "full" CDF,
#   the above approach would all of a sudden return a different full name!
# o Added clearCache().
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
