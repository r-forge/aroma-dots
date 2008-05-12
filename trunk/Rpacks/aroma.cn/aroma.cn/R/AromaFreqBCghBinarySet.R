###########################################################################/**
# @RdocClass AromaFreqBCghBinarySet
#
# @title "The AromaFreqBCghBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaFreqBCghBinarySet object represents a set of 
#  @see "AromaFreqBCghBinaryFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaFreqBCghBinarySet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaFreqBCghBinarySet");
})


setMethodS3("byName", "AromaFreqBCghBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cghData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    fromFiles(static, path=path, ..., pattern=".*,freqB[.]asb$");
  })
}, static=TRUE) 




############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
