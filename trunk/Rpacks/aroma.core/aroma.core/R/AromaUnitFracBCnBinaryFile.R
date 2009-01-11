###########################################################################/**
# @RdocClass AromaUnitFracBCnBinaryFile
#
# @title "The AromaUnitFracBCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitFracBCnBinaryFile is a @see "AromaUnitTabularBinaryFile".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaUnitFracBCnBinaryFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), "AromaUnitFracBCnBinaryFile"
  );
})



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
