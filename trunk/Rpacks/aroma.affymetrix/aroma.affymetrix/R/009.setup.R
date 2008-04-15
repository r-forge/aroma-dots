.setupAromaAffymetrix <- function(pkg, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch base::matrix() to be more memory efficient when 'dimnames==NULL'.
  .patchMatrix();

  # Patch log2()/log10() that are slow to display warnings
  .patchLog2();

  # Patch affxparser::findCdf()
  .patchCdfMergeStrands();
#  reassignInPackage("findFiles", "affxparser", findFiles.patch);
#  reassignInPackage("findCdf", "affxparser", findCdf.patch);



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

  template <- list(
    memory = list(
      gcArrayFrequency = 50
    ),

    rules = list(
      allowAsciiCdfs = FALSE
    ),

    output = list(
      # Max number of arrays for which to report timestamps
      timestampsThreshold = 500
    ),

    models = list(
      RmaPlm = list( 
       # Number of cells *and* arrays for using median polish
        medianPolishThreshold  = c( 500, 6),
       # Number of cells *and* arrays for skipping unit group
        skipThreshold          = c(5000, 1)
      )
    ),

    annotationData = list(
      paths = list()
    ),

    system = list(
      checkForUpdates = FALSE,
      checkInterval = "onEachLoad"
    )
  );

  # Copy settings from template, if missing
  for (dir in names(template)) {
    if (!dir %in% names(settings))
      settings[[dir]] <- list();
    for (subdir in names(template[[dir]])) {
      if (!subdir %in% names(settings[[dir]]))
        settings[[dir]][[subdir]] <- template[[dir]][[subdir]];
    }
  }

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
  # Check for updates and patches?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  settings <- getOption("aroma.affymetrix.settings");
  if (identical(settings$system$checkForUpdates, TRUE)) {
    interval <- settings$system$checkInterval;
    if (is.null(interval))
      interval <- "onEachLoad";

    if (identical(interval, "onEachLoad")) {
      if (interactive()) {
        update(pkg, verbose=TRUE);
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that digest() gives a consistent result across R versions
  # and platforms.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!identical(getOption("aroma.affymetrix::assertDigest"), FALSE)) {
    .assertDigest("error");
  }
} # .setupAromaAffymetrix()



############################################################################
# HISTORY:
# 2008-02-14
# o Renamed existing threshold hold to 'timestampsThreshold', 
#   'medianPolishThreshold', and 'skipThreshold'.
# 2008-02-12
# o Added default values for settings 'models$RmaPlm$...'.
# 2008-01-30
# o Added default values for settings 'rules$allowAsciiCdfs' and
#   'output$maxNbrOfArraysForTimestamps'.
# 2007-12-13
# o Added code for automatic updates on startup.  In active by default.
# o Added settings for 'checkForPatches' and 'checkInterval'.
# o Now the settings are set according to a tempate, if missing.
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
