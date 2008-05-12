###########################################################################/**
# @RdocClass AromaFreqBCnBinarySet
#
# @title "The AromaFreqBCnBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaFreqBCnBinarySet object represents a set of 
#  @see "AromaFreqBCnBinaryFile"s with \emph{identical} chip types.
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
setConstructorS3("AromaFreqBCnBinarySet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaFreqBCnBinarySet");
})


setMethodS3("byName", "AromaFreqBCnBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cnData") {
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
