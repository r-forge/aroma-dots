.setupAromaAffymetrix <- function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch base::matrix() to be more memory efficient when 'dimnames==NULL'.
  .patchMatrix();

  # Package as.<basic data type>()? (for memory efficiencies)
  .patchAsDataTypes();

  # Make sure there is rowMedians() supporting missing values
  .patchRowMedians();

  # Patch affxparser::findCdf()
  reassignInPackage("findCdf", "affxparser", findCdf.patch);



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




  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that digest() gives a consistent result across R versions
  # and platforms.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!identical(getOption("aroma.affymetrix::assertDigest"), FALSE))
    .assertDigest("error");
} # .setupAromaAffymetrix()


############################################################################
# HISTORY:
# 2007-08-30
# o Added "patch" to make sure that there is rowMedians() supporting 
#   missing values.
# 2007-07-04
# o Removed the patch for digest(); digest v0.3.0 solved the problem.
# o Added a patch of functions in 'base', e.g. matrix().
# 2007-04-04
# o Moved the patch of digest() here.
# 2007-03-07
# o Added test for consistency of digest().
# 2007-03-06
# o Added onLoad hook function for library() and require() to call
#   fixSearchPath() of the package, which reorders the search path so that
#   problematic packages are after this package in the search path.
# 2007-02-22
# o Added default settings stubs.
# 2007-02-12
# o Created.
############################################################################
