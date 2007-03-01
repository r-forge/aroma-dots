.setupAromaAffymetrix <- function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch some of the affxparser function (for now)
  .patchAffxparser();

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Customize affxparser
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add custom findCdf() function to affxparser.  This is need to be
  # able to locate CDFs in annotationData/chipTypes/<chipType>/.
  setCustomFindCdf(function(...) {
    AffymetrixCdfFile$findByChipType(..., .useAffxparser=FALSE);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Package settings (settings might change)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load
  settings <- getOption("aroma.affymetrix.settings");	
  if (is.null(settings))
    settings <- list();

  if (!"annotationData" %in% names(settings))
    settings$annotationData <- list();
  
#  if (!"aliases" %in% names(settings))
#    settings$annotationData$aliases <- list();
  
  if (!"paths" %in% names(settings))
    settings$annotationData$paths <- list();

  options("aroma.affymetrix.settings"=settings);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply downloaded patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  patchPackage("aroma.affymetrix");
} # .setupAromaAffymetrix()

############################################################################
# HISTORY:
# 2007-02-22
# o Added default settings stubs.
# 2007-02-12
# o Created.
############################################################################
