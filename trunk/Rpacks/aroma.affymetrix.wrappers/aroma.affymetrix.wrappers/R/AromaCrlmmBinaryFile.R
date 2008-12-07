###########################################################################/**
# @RdocClass AromaCrlmmBinaryFile
#
# @title "The AromaCrlmmBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaCrlmmBinaryFile is a @see "aroma.core::AromaTabularBinaryFile".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.core::AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaCrlmmBinaryFile", function(...) {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  extend(AromaSignalBinaryFile(...), "AromaCrlmmBinaryFile"
  );
})


setMethodS3("allocate", "AromaCrlmmBinaryFile", function(static, ..., nbrOfStrands=2, types=c("integer", rep("double", 1+3*nbrOfStrands)), sizes=c(1, rep(4, 1+3*nbrOfStrands))) { 
  allocate.AromaSignalBinaryFile(static, types=types, sizes=sizes, ...);
})



############################################################################
# HISTORY:
# 2008-12-05
# o Created.
############################################################################
