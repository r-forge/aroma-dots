.setupAromaAffymetrix <- function(...) {
  # Patch some of the affxparser function (for now)
  .patchAffxparser();
  
  # Add custom findCdf() function to affxparser.  This is need to be
  # able to locate CDFs in annotationData/chipTypes/<chipType>/.
  setCustomFindCdf(function(...) {
    AffymetrixCdfFile$findByChipType(..., .useAffxparser=FALSE);
  });
} # .setupAromaAffymetrix()

############################################################################
# HISTORY:
# 2007-02-12
# o Created.
############################################################################
