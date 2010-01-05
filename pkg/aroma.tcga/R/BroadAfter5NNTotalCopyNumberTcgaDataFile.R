setConstructorS3("BroadAfter5NNTotalCopyNumberTcgaDataFile", function(...) {
  extend(BroadTotalCopyNumberTcgaDataFile(...), "BroadAfter5NNTotalCopyNumberTcgaDataFile");
})


setMethodS3("getExtensionPattern", "BroadAfter5NNTotalCopyNumberTcgaDataFile", function(this, default="[.](after_5NN[.]copynumber[.]data[.]txt)$", ...) {
  NextMethod("getExtensionPattern", this, default=default, ...);
}, static=TRUE)



############################################################################
# HISTORY:
# 2010-01-04
# o Created from BroadTotalCopyNumberTcgaDataFile.R.
############################################################################
