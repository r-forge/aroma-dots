###########################################################################/**
# @RdocClass RawMultisourceCopyNumberModel
#
# @title "The RawMultisourceCopyNumberModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents an idenity multi-source copy-number model which
#  returns the input as is.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Passed to the constructor of the superclass.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RawMultisourceCopyNumberModel", function(...) {
  extend(MultisourceCopyNumberChromosomalModel(...), "RawMultisourceCopyNumberModel");
})

setMethodS3("getSetTag", "RawMultisourceCopyNumberModel", function(this, ...) {
  "raw";
})



##############################################################################
# HISTORY:
# 2009-01-25
# o Created from RawCopyNumberModel.
##############################################################################
