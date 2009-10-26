setConstructorS3("TcgaDataFile", function(...) {
  extend(MageTabDataMatrixFile(...), "TcgaDataFile");
})

setMethodS3("getExtensionPattern", "TcgaDataFile", function(static, ...) {
  "[.](data[.]txt)$";
}, static=TRUE)



############################################################################
# HISTORY:
# 2009-08-23
# o Created.
############################################################################
