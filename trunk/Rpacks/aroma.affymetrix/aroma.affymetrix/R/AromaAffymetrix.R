setConstructorS3("AromaAffymetrix", function(...) {
  extend(Package("aroma.affymetrix"), "AromaAffymetrix");
})

setMethodS3("fixSearchPath", "AromaAffymetrix", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Making sure the search path is compatible with ", getName(this));

  toPath <- sprintf("package:%s", getName(this));

  verbose && cat(verbose, "Search path before:");
  verbose && print(verbose, search());

  # Problematic package that must be after this package on the search path
  pkgs <- c("affy", "affyPLM", "EBImage");

  # Move those package, if they are loaded.
  pkgsMoved <- c();
  for (pkg in pkgs) {
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

############################################################################
# HISTORY:
# 2007-03-06
# o Created.
############################################################################
