###########################################################################/**
# @RdocClass AromaTotalCghBinaryFile
#
# @title "The AromaTotalCghBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaTotalCghBinaryFile is a @see "AromaTabularBinaryFile".
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
setConstructorS3("AromaTotalCghBinaryFile", function(...) {
  extend(AromaSignalBinaryFile(...), "AromaTotalCghBinaryFile"
  );
})



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
