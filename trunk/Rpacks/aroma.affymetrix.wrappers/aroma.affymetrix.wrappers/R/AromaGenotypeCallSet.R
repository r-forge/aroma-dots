###########################################################################/**
# @RdocClass AromaGenotypeCallSet
#
# @title "The AromaGenotypeCallSet class"
#
# \description{
#  @classhierarchy
#
#  An AromaGenotypeCallSet object represents a set of 
#  @see "AromaCrlmmBinaryFile"s with \emph{identical} chip types.
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
setConstructorS3("AromaGenotypeCallSet", function(...) {
  extend(AromaSignalBinarySet(...), "AromaGenotypeCallSet");
})


setMethodS3("byName", "AromaGenotypeCallSet", function(static, name, tags=NULL, ..., chipType=NULL, paths="crlmmData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  fromFiles(static, path=path, ...);
}, static=TRUE) 


setMethodS3("fromFiles", "AromaGenotypeCallSet", function(static, ...) {
  suppressWarnings({
    fromFiles.GenericDataFileSet(static, ..., pattern=".*[.]agc$");
  })
})

setMethodS3("findUnitsTodo", "AromaGenotypeCallSet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  ce <- getFile(this, length(this));
  findUnitsTodo(ce, ...);
})


setMethodS3("extractCalls", "AromaGenotypeCallSet", function(this, units=NULL, ...) {
  values <- extractMatrix(this, rows=units, column=1, ...);
  isNA <- whichVector(values == as.integer(255));
  values[isNA] <- as.integer(NA);
  values;
})


setMethodS3("extractConfidenceScores", "AromaGenotypeCallSet", function(this, units=NULL, ...) {
  values <- extractMatrix(this, rows=units, column=2, ...);
  values;
})

############################################################################
# HISTORY:
# 2008-12-09
# o Created.
############################################################################
