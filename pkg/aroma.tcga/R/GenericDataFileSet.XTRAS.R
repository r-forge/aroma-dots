setMethodS3("appendFullNameTranslatorBySdrfFile", "FullNameInterface", function(this, sdrf, ...) {
  sdrfs <- SdrfFileSet(sdrf);
  appendFullNameTranslator(this, sdrfs);
})


setMethodS3("appendFullNameTranslatorBySdrfFileSet", "FullNameInterface", function(this, set, ...) {
  # Argument 'set':
  className <- "SdrfFileSet";
  if (!inherits(set, className)) {
    throw("Argument 'sdrf' is not of class ", className, ": ", class(set)[1]);
  }

  # Build a fullnames translator function based on the SDRF
  fcn <- makeFullNamesTranslator(set, ...);

  # Apply the translator function
  appendFullNameTranslator(this, fcn);
})


############################################################################
# HISTORY:
# 2009-10-23
# o OPTIMIZATION: Now makeFullNamesTranslator() of SdrfFileSet caches the
#   result in the memory so it does not have to reread the files each times.
# 2009-10-02
# o Added setFullNamesTranslatorBySdrfFileSet() to GenericDataFileSet.
# 2009-10-01
# o Added setFullNamesTranslatorBySdrfFile() to GenericDataFileSet.
# o Created.
############################################################################
