###########################################################################/**
# @RdocClass AromaCrlmmBinarySet
#
# @title "The AromaCrlmmBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaCrlmmBinarySet object represents a set of 
#  @see "AromaCrlmmBinaryFile"s with \emph{identical} chip types.
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
setConstructorS3("AromaCrlmmBinarySet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaCrlmmBinarySet");
})


setMethodS3("byName", "AromaCrlmmBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cnData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    fromFiles(static, path=path, ..., pattern=".*,CRLMM[.]acu$");
  })
}, static=TRUE) 




############################################################################
# HISTORY:
# 2008-12-05
# o Created.
############################################################################
