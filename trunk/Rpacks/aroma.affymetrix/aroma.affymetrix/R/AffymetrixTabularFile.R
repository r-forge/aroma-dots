setConstructorS3("AffymetrixTabularFile", function(...) {
  extend(GenericTabularFile(...), "AffymetrixTabularFile");
})



setMethodS3("findByChipType", "AffymetrixTabularFile", function(static, chipType, tags=NULL, pattern=NULL, ...) {
  if (is.null(pattern)) {
    name <- paste(c(chipType, tags), collapse=",");
    pattern <- sprintf("^%s.*[.]...$", name);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern, ...);
  pathname;
}, static=TRUE, protected=TRUE)



setMethodS3("fromChipType", "AffymetrixTabularFile", function(static, chipType, ...) {
  # Search for the genome information file
  pathname <- static$findByChipType(chipType, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix tabular file: ", chipType);
  newInstance(static, pathname);
}, static=TRUE)



############################################################################
# HISTORY:
# 2007-09-14
# o Now inheriting from GenericTabularFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
