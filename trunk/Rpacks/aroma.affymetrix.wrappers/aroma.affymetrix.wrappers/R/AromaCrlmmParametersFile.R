###########################################################################/**
# @RdocClass AromaCrlmmParametersFile
#
# @title "The AromaCrlmmParametersFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaCrlmmParametersFile is a @see "aroma.core::AromaTabularBinaryFile".
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
setConstructorS3("AromaCrlmmParametersFile", function(...) {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  extend(AromaSignalBinaryFile(...), "AromaCrlmmParametersFile"
  );
})


setMethodS3("allocate", "AromaCrlmmParametersFile", function(static, ..., nbrOfStrands=2, types=rep("double", 1+3*nbrOfStrands), sizes=rep(4, 1+3*nbrOfStrands), signed=rep(TRUE, 1+3*nbrOfStrands)) { 
  res <- allocate.AromaSignalBinaryFile(static, types=types, sizes=sizes, signed=signed, ...);
  res;
})



setMethodS3("findUnitsTodo", "AromaCrlmmParametersFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Identifying non-fitted units in file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));

  # Reading all calls
  values <- this[,1,drop=TRUE];

  units <- whichVector(values == 0);
  verbose && exit(verbose);

  units;
})




############################################################################
# HISTORY:
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-05
# o Created.
############################################################################
