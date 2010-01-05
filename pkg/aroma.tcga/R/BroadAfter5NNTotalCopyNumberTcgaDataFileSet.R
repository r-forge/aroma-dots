setConstructorS3("BroadAfter5NNTotalCopyNumberTcgaDataFileSet", function(...) {
  extend(BroadTotalCopyNumberTcgaDataFileSet(...), "BroadAfter5NNTotalCopyNumberTcgaDataFileSet");
})


setMethodS3("exportTotalAndFracB", "BroadAfter5NNTotalCopyNumberTcgaDataFileSet", function(this, tags=c("*", "5NN"), ...) {
  NextMethod("exportTotalAndFracB", this, tags=tags, ...);
})


############################################################################
# HISTORY:
# 2010-01-04
# o Created.
############################################################################
