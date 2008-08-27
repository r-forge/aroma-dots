setConstructorS3("AromaAffymetrix", function(...) {
  extend(Package("aroma.affymetrix"), "AromaAffymetrix");
})

setMethodS3("fixSearchPath", "AromaAffymetrix", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Making sure the search path is compatible with ", getName(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RULES
  # 2008-02-27:
  # o affy must be after R.huge, otherwise the former overrides the
  #   generic colnames() function of the latter.
  # o affyPLM must be after aroma.affymetrix.
  # o EBImage must be after aroma.affymetrix.
  # 2008-08-27:
  # o affy must be after aroma.light, otherwise the former overrides
  #   the generic plotDensity() function of the latter.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Figure out which of aroma.affymetrix, R.huge, and aroma.light is
  # last on the search path. (Since aroma.affymetrix loads both R.huge
  # and aroma.light, it is guaranteed to be on the search path before
  # those, but to avoid bugs in the future we don't assume that).
  pkgs <- c("aroma.affymetrix", "aroma.light", "R.huge");
  idxs <- match(sprintf("package:%s", pkgs), search());
  lastPkg <- pkgs[which.max(idxs)];
  toPath <- sprintf("package:%s", lastPkg);

  verbose && cat(verbose, "Search path before:");
  verbose && print(verbose, search());

  # Problematic package that must be after this package on the search path
  pkgsToMove <- c("affy", "affyPLM", "EBImage");

  # Move those package, if they are loaded.
  pkgsMoved <- c();
  for (pkg in pkgsToMove) {
    path <- sprintf("package:%s", pkg);
    if (path %in% search()) {
      # Need to move?
      from <- match(path, search());
      to <- match(toPath, search());
      if (from < to) {
        verbose && printf(verbose, "Moving package: %s (%s)\n", pkg, path);
        pkgsMoved <- c(pkgsMoved, 
                    moveInSearchPath(from=path, to=toPath, where="after"));
      }
    }
  }

  verbose && cat(verbose, "Search path after:");
  verbose && print(verbose, search());

  verbose && exit(verbose);

  # Return the package actually moved
  invisible(pkgsMoved);
})


setMethodS3("update", "AromaAffymetrix", function(object, patch=TRUE, ..., verbose=FALSE) {
  # To please R CMD check
  this <- object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Checking for and install updates");
  verbose && cat(verbose, "Package: ", getName(this));
  verbose && printf(verbose, "Current version: v%s (%s)\n", 
                                        getVersion(this), getDate(this));


  state <- 0;

  url <- "http://www.braju.com/R/hbLite.R";
  verbose && enter(verbose, "Trying to download update script");
  verbose && cat(verbose, "URL: ", url);
  hbInstall <- NULL;
  tryCatch({
    suppressWarnings({
      source(url);
    })
    state <- 1;
    verbose && exit(verbose);
  }, error = function(ex) {
    verbose && exit(verbose, suffix="failed");
    throw(ex);
  })

  verbose && enter(verbose, "Launching update command");
  verbose && printf(verbose, "Call: hbInstall(\"%s\")", getName(this));
  tryCatch({
    hbInstall(getName(this));
    state <- 2;
  }, error = function(ex) {
    verbose && exit(verbose, suffix="failed");
    throw(ex);
  })

  verbose && cat(verbose, "Package has been updated.");

  if (patch) {
    patch(this, ..., verbose=verbose);
  }

  verbose && exit(verbose);

  invisible(state);
})



setMethodS3("patch", "AromaAffymetrix", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Checking for and install patches");
  verbose && cat(verbose, "Package: ", getName(this));
  verbose && printf(verbose, "Current version: v%s (%s)\n", 
                                        getVersion(this), getDate(this));


  state <- 0;

  verbose && enter(verbose, "Trying to download patches");
  tryCatch({
    downloadPackagePatch(getName(this), verbose=verbose);
    state <- 1;
    verbose && exit(verbose);
  }, error = function(ex) {
    verbose && exit(verbose, suffix="failed");
    throw(ex);
  })

  verbose && cat(verbose, "Package has been patched.");
  
  verbose && exit(verbose);

  invisible(state);
})


############################################################################
# HISTORY:
# 2008-08-27
# o Now the affy, affyPLM, and EBImage package are forced to be after all
#   of aroma.affymetrix, aroma.light, and R.huge on the search() path.
# 2007-12-13
# o Added update() and patch() to the AromaAffymetrix Package class.
# 2007-03-06
# o Created.
############################################################################
