setConstructorS3("AromaAffymetrixWrappers", function(...) {
  extend(Package("aroma.affymetrix.wrappers"), "AromaAffymetrixWrappers");
})

setMethodS3("fixSearchPath", "AromaAffymetrixWrappers", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Making sure the search path is compatible with ", getName(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RULES
  # 2008-02-27:
  # o oligo must be after aroma.affymetrix.wrappers.
  # o affyPLM must be after aroma.affymetrix.wrappers.
  # o EBImage must be after aroma.affymetrix.wrappers.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  toPath <- sprintf("package:%s", getName(this));

  verbose && cat(verbose, "Search path before:");
  verbose && print(verbose, search());

  # Problematic package that must be after this package on the search path
  pkgsToMove <- c("affy", "affyPLM", "oligo");
  pkgsToMove <- c();

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


############################################################################
# HISTORY:
# 2008-12-05
# o Created from aroma.affymetrix.
############################################################################
