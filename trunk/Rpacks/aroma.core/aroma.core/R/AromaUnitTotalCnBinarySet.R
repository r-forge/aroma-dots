###########################################################################/**
# @RdocClass AromaUnitTotalCnBinarySet
#
# @title "The AromaUnitTotalCnBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitTotalCnBinarySet object represents a set of 
#  @see "AromaUnitTotalCnBinaryFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaUnitTotalCnBinarySet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "AromaUnitTotalCnBinarySet");
})


setMethodS3("byName", "AromaUnitTotalCnBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cnData") {
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
