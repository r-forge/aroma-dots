setMethodS3("appendFullNameTranslatorBySdrfFile", "FullNameInterface", function(this, sdrf, ...) {
  sdrfs <- SdrfFileSet(sdrf);
  appendFullNameTranslator(this, sdrfs);
})


setMethodS3("appendFullNameTranslatorBySdrfFileSet", "FullNameInterface", function(this, set, ...) {
  # Argument 'set':
  set <- Arguments$getInstanceOf(set, "SdrfFileSet");

  # Build a fullnames translator function based on the SDRF
  fcn <- makeFullNamesTranslator(set, ...);

  # Apply the translator function
  appendFullNameTranslator(this, fcn);
})


setMethodS3("extractByTcgaType", "GenericDataFileSet", function(this, pattern, ...) {
  getTypeTag <- function(this, ...) {
    gsub("[A-Z]$", "", getTags(this)[1]);
  }

  # Argument 'typePattern':
  pattern <- Arguments$getRegularExpression(pattern);

  # Sanity check
  names <- getFullNames(this);
  patterns <- BiospecimenCoreResource$getBarcodePatterns();
  pattern2 <- sprintf("^%s,%s", patterns$patient, patterns$sampleId);
  nok <- (regexpr(pattern2, names) == -1);
  if (any(nok)) {
    throw(sprintf("%s '%s' contains %s files with non-TCGA names.", 
          class(this)[1], getFullName(this), sum(nok)));
  }

  # Identify sample types matching the pattern
  types <- sapply(this, getTypeTag);
  keep <- grep(pattern, types);

  # Extract matching samples
  res <- extract(this, keep);

  res;
})


############################################################################
# HISTORY:
# 2010-01-05 
# o Added extractByTcgaType().
# 2009-10-23
# o OPTIMIZATION: Now makeFullNamesTranslator() of SdrfFileSet caches the
#   result in the memory so it does not have to reread the files each times.
# 2009-10-02
# o Added setFullNamesTranslatorBySdrfFileSet() to GenericDataFileSet.
# 2009-10-01
# o Added setFullNamesTranslatorBySdrfFile() to GenericDataFileSet.
# o Created.
############################################################################
