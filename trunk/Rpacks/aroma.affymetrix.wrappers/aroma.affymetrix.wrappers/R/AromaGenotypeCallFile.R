###########################################################################/**
# @RdocClass AromaGenotypeCallFile
#
# @title "The AromaGenotypeCallFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaGenotypeCallFile is a @see "aroma.core::AromaTabularBinaryFile".
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
setConstructorS3("AromaGenotypeCallFile", function(...) {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  extend(AromaSignalBinaryFile(...), "AromaGenotypeCallFile"
  );
})


setMethodS3("allocate", "AromaGenotypeCallFile", function(static, ..., types=c("integer", "double"), sizes=c(1, 4), signed=c(FALSE, TRUE)) { 
  res <- allocate.AromaSignalBinaryFile(static, types=types, sizes=sizes, signed=signed, ...);
  # Default call is NA (=255)
  res[,1] <- as.integer(255);
  res;
})



setMethodS3("findUnitsTodo", "AromaGenotypeCallFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
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


setMethodS3("extractCalls", "AromaGenotypeCallFile", function(this, units=NULL, ..., verbose=FALSE) {
  values <- extractMatrix(this, rows=units, column=1, ...);
  isNA <- whichVector(values == as.integer(255));
  values[isNA] <- as.integer(NA);
  values;
})


setMethodS3("extractConfidenceScores", "AromaGenotypeCallFile", function(this, units=NULL, ..., verbose=FALSE) {
  values <- extractMatrix(this, rows=units, column=2, ...);
  values;
})


############################################################################
# HISTORY:
# 2008-12-08
# o Recreated.  Now only with columns (genotypeCall, confidenceScore).
# o Added findUnitsTodo() and extractCalls().
# 2008-12-05
# o Created.
############################################################################
