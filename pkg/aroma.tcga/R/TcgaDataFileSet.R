setConstructorS3("TcgaDataFileSet", function(...) {
  extend(MageTabFileSet(...), "TcgaDataFileSet");
})

setMethodS3("getDefaultFullName", "TcgaDataFileSet", function(this, parent=1, ...) {
  NextMethod("getDefaultFullName", this, parent=1, ...);
})


############################################################################
# HISTORY:
# 2010-01-04
# o Added getDefaultFullName(... parent=1) to TcgaDataFileSet so that
#   it handles the updated ditto in R.filesets (now with parent=0).
# 2009-10-25
# o Created.
############################################################################
