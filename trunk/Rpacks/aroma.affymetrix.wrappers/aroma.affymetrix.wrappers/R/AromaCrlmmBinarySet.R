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



setMethodS3("findUnitsTodo", "AromaCrlmmBinarySet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  ce <- getFile(this, length(this));
  findUnitsTodo(ce, ...);
})


setMethodS3("extractCalls", "AromaCrlmmBinarySet", function(this, units=NULL, ...) {
  values <- extractMatrix(this, rows=units, column=1, ...);
  isNA <- whichVector(values == as.integer(255));
  values[isNA] <- as.integer(NA);
  values;
})

############################################################################
# HISTORY:
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-06
# o Created.
############################################################################
