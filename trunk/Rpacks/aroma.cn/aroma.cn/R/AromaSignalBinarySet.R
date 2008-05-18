###########################################################################/**
# @RdocClass AromaSignalBinarySet
#
# @title "The AromaSignalBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaSignalBinarySet object represents a set of 
#  @see "AromaSignalBinaryFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaSignalBinarySet", function(...) {
  extend(AromaTabularBinarySet(...), "AromaSignalBinarySet");
})


setMethodS3("findByName", "AromaSignalBinarySet", function(static, ..., chipType=NULL) {
  # NextMethod() does not work here.
  findByName.GenericDataFileSet(static, ..., subdirs=chipType);
}, static=TRUE) 


setMethodS3("byName", "AromaSignalBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths=NULL, pattern="[.]asb$") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    fromFiles(static, path=path, ..., pattern=pattern);
  })
}, static=TRUE) 



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("validate", "AromaSignalBinarySet", function(this, ...) {
  chipTypes <- lapply(this, FUN=getChipType);
  chipTypes <- unique(chipTypes);
  if (length(chipTypes) > 1) {
    throw("The located ", class(this)[1], " contains files with different chip types: ", paste(chipTypes, collapse=", "));
  }

  NextMethod("validate", this, ...);
}, protected=TRUE)


setMethodS3("getPlatform", "AromaSignalBinarySet", function(this, ...) {
  getPlatform(getFile(this, 1), ...);
})


setMethodS3("getChipType", "AromaSignalBinarySet", function(this, ...) {
  getChipType(getFile(this, 1), ...);
}) 


setMethodS3("getAromaUgpFile", "AromaSignalBinarySet", function(this, ...) {
  getAromaUgpFile(getFile(this,1), ...);
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



############################################################################
# HISTORY:
# 2008-05-11
# o Created.
############################################################################
