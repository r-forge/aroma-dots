setMethodS3("cleanToRealChipType", "AffymetrixCdfFile", function(static, cleanCdfName, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # We have to turn to an ad hoc solution.

  #  1) Find all unique known CDF files
  cdfFiles <- static$findByChipType(firstOnly=FALSE);
  cdfNames <- basename(cdfFiles);
  cdfNames <- gsub("[.](c|C)(d|D)(f|F)$", "", cdfNames);
  cdfNames <- gsub("[.]text$", "", cdfNames);
  cdfNames <- unique(cdfNames);

  #  2) Get the clean names of these
  cleanNames <- sapply(cdfNames, FUN=affy::cleancdfname, addcdf=FALSE);

  #  3) Find the one match the one we are looking for
  pattern <- sprintf("^%s$", cleanCdfName);
  cdfNames <- cdfNames[grep(pattern, cleanNames)];
  if (length(cdfNames) == 0)
    cdfNames <- NULL;
  cdfNames;
}, static=TRUE, private=TRUE)



setMethodS3("fromCleanChipType", "AffymetrixCdfFile", function(static, cleanChipType, ...) {
  chipType <- cleanToRealChipType(static, cleanChipType, ...);
  fromChipType(static, chipType);
}, static=TRUE, private=TRUE)


setMethodS3("getOligoToCelMap", "AffymetrixCdfFile", function(this, invertMap=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving the feature-to-CEL index map");
  verbose && cat(verbose, "Chip type: ", getChipType(this));

  verbose && enter(verbose, "Retrieving (X,Y) from the platform-design package");
  pd <- PlatformDesign(this);
  x <- getFeatureInfo(pd, "X")-1;
  y <- getFeatureInfo(pd, "Y")-1;
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculate the map");
  map <- as.integer(nbrOfColumns(this)*y+x+1);
  rm(x, y);
  if (!invertMap) {
    # This negation is correct. /HB 2006-12-07
    map <- affxparser::invertMap(map);
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  map;
}, private=TRUE)


setMethodS3("getPlatformDesignObject", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  pd <- PlatformDesign(this);

  verbose && enter(verbose, "Retrieves the platform-design object");
  pdo <- getPlatformDesignObject(pd);
  verbose && exit(verbose);

  attr(pd, "name") <- getPackageName(pd);
  attr(pd, "chipType") <- getChipType(pd);
  attr(pd, "cleanChipType") <- getChipType(pd, clean=TRUE);

  pd;
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-06-11
# o BUG FIX: Called getPlatformDesign() instead of getPlatformDesignObject()
#   for a PlatformDesign object in getPlatformDesignObject() of 
#   AffymetrixCdfFile.
# 2006-12-07
# o Added getOligoToCelMap() and getPlatformDesign().
# 2006-12-06
# o Created.
############################################################################
