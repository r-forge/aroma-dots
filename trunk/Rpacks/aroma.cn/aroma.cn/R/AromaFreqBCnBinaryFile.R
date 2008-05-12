###########################################################################/**
# @RdocClass AromaFreqBCnBinaryFile
#
# @title "The AromaFreqBCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaFreqBCnBinaryFile is a @see "AromaTabularBinaryFile".
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
setConstructorS3("AromaFreqBCnBinaryFile", function(...) {
  extend(AromaSignalBinaryFile(...), "AromaFreqBCnBinaryFile"
  );
})



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
