###########################################################################/**
# @RdocClass AromaTotalCnBinaryFile
#
# @title "The AromaTotalCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaTotalCnBinaryFile is a @see "AromaTabularBinaryFile".
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
setConstructorS3("AromaTotalCnBinaryFile", function(...) {
  extend(AromaSignalBinaryFile(...), "AromaTotalCnBinaryFile"
  );
})



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
