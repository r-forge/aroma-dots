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


setMethodS3("allocate", "AromaCrlmmBinaryFile", function(static, ..., nbrOfStrands=2, types=c("integer", rep("double", 1+3*nbrOfStrands)), sizes=c(1, rep(4, 1+3*nbrOfStrands)), signed=c(FALSE, rep(TRUE, 1+3*nbrOfStrands))) { 
  res <- allocate.AromaSignalBinaryFile(static, types=types, sizes=sizes, signed=signed, ...);
  # Default call is NA (=255)
  res[,1] <- as.integer(255);
  res;
})



setMethodS3("findUnitsTodo", "AromaCrlmmBinaryFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
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
  calls <- this[,1,drop=TRUE];

  # Locate the non-fitted ones
  isNA <- (calls == as.integer(255));

  units <- whichVector(isNA);
  rm(isNA);
  verbose && exit(verbose);

  units;
})


setMethodS3("extractCalls", "AromaCrlmmBinaryFile", function(this, units=NULL, ..., verbose=FALSE) {
  values <- extractMatrix(this, rows=units, column=1, ...);
  isNA <- whichVector(values == as.integer(255));
  values[isNA] <- as.integer(NA);
  values;
})


############################################################################
# HISTORY:
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-05
# o Created.
############################################################################
