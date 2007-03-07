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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fix the search path every time a package is loaded
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setHook("base::library:onLoad", function(...) {
    # Fix the search path
    pkgs <- fixSearchPath(aroma.affymetrix);
    if (length(pkgs) > 0) {
      warning("Packages reordered in search path: ", 
                                            paste(pkgs, collapse=", "));
    }
  }, action="append");
} # .setupAromaAffymetrix()


############################################################################
# HISTORY:
# 2007-03-06
# o Added onLoad hook function for library() and require() to call
#   fixSearchPath() of the package, which reorders the search path so that
#   problematic packages are after this package in the search path.
# 2007-02-22
# o Added default settings stubs.
# 2007-02-12
# o Created.
############################################################################
