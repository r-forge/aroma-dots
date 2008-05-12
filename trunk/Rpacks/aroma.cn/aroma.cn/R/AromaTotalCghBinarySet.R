###########################################################################/**
# @RdocClass AromaTotalCghBinarySet
#
# @title "The AromaTotalCghBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaTotalCghBinarySet object represents a set of 
#  @see "AromaTotalCghBinaryFile"s with \emph{identical} chip types.
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
setConstructorS3("AromaTotalCghBinarySet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaTotalCghBinarySet");
})


setMethodS3("byName", "AromaTotalCghBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cghData") {
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
