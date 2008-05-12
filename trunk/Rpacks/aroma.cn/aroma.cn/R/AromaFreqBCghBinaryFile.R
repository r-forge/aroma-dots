###########################################################################/**
# @RdocClass AromaFreqBCghBinaryFile
#
# @title "The AromaFreqBCghBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaFreqBCghBinaryFile is a @see "AromaTabularBinaryFile".
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
#*/########################################################################### 
setConstructorS3("AromaFreqBCghBinaryFile", function(...) {
  extend(AromaSignalBinaryFile(...), "AromaFreqBCghBinaryFile"
  );
})



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
