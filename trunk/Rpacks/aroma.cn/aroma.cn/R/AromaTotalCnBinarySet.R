###########################################################################/**
# @RdocClass AromaTotalCnBinarySet
#
# @title "The AromaTotalCnBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaTotalCnBinarySet object represents a set of 
#  @see "AromaTotalCnBinaryFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.core::AromaSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaTotalCnBinarySet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaTotalCnBinarySet");
})


setMethodS3("byName", "AromaTotalCnBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cnData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    fromFiles(static, path=path, ..., pattern=".*,total[.]asb$");
  })
}, static=TRUE) 




############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
