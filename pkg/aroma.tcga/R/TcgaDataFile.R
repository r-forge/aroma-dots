setConstructorS3("TcgaDataFile", function(...) {
  extend(MageTabDataMatrixFile(...), "TcgaDataFile");
})

setMethodS3("getExtensionPattern", "TcgaDataFile", function(this, default="[.](data[.]txt)$", ...) {
  NextMethod("getExtensionPattern", this, default=default, ...);  
}, static=TRUE)



############################################################################
# HISTORY:
# 2010-01-04
# o Made the default of getExtensionPattern() for TcgaDataFile an argument.
# 2009-08-23
# o Created.
############################################################################
