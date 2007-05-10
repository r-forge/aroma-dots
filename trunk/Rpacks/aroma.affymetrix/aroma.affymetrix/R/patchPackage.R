###########################################################################/**
# @RdocDefault patchPackage
#
# @title "Applies patches for a specific package"
#
# \description{
#  @get "title" from the package online reprocitory.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{pkgName}{The name of the package to be patched."}
#   \item{paths}{A @character @vector of paths containing package patches.}
#   \item{deleteOld}{If @TRUE, old patch directories are deleted.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns the number of files sourced.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("patchPackage", "default", function(pkgName, paths="patches", deleteOld=TRUE, verbose=FALSE, ...) {
  require(R.utils) || stop("Package not loaded: R.utils");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  findPatchDirectories <- function(pkgName, path, verbose=FALSE, ...) {
    # Locate patch root directory
    rootPath <- file.path(path, pkgName);
    if (!isDirectory(rootPath)) {
      verbose && cat(verbose, "No patch root directory found: ", rootPath);
      return(c());
    }

    verbose && enter(verbose, "Root path: ", rootPath);
  
    # Search for patch directories
    pattern <- "^20[0-9][0-9][01][0-9][0-3][0-9](|[a-z])";
    paths <- list.files(path=rootPath, pattern=pattern, full.names=TRUE);
    if (length(paths) > 0)
      paths <- paths[sapply(paths, FUN=isDirectory)];
    if (length(paths) == 0) {
      verbose && cat(verbose, "No patch directory found in root path: ", rootPath);
      return(c());
    }
  
    paths;
  } # findPatchDirectories()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'deleteOld':
  deleteOld <- Arguments$getVerbose(deleteOld);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Patching package ", pkgName);

#  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);


  # 1. Scan for patch directories
  paths <- lapply(paths, FUN=function(path) {
    findPatchDirectories(pkgName, path=path, verbose=verbose);
  })
  paths <- unlist(paths, use.names=FALSE);

  # Nothing do to?
  if (length(paths) == 0) {
    verbose && exit(verbose);
    return(invisible(0));
  }

  # 2. Order paths by date
  o <- order(basename(paths), decreasing=TRUE);
  paths <- paths[o];
  verbose && cat(verbose, "Identified patch paths:");
  verbose && print(verbose, paths);

  # Get package date
  pkgDate <- packageDescription("aroma.affymetrix")$Date;
  pkgDate <- gsub("-", "", pkgDate);
  verbose && cat(verbose, "Package date: ", pkgDate);

  # 3. Keep only the most recent one, if it is newer than the package
  path <- paths[1];
  if (basename(path) < pkgDate)
    path <- NULL;

  # 4. Remove all older patches?
  if (deleteOld) {
    # Delete patch paths that are older than the package version
    delPaths <- paths[basename(paths) < pkgDate];

    # Delete also patch paths that are older than the latest patch
    delPaths <- c(delPaths, paths[-1]);

    # Any paths to be removed?
    delPaths <- unique(delPaths);
    if (length(delPaths) > 0) {
      verbose && enter(verbose, "Deleting deprecated patch directories");
      res <- sapply(delPaths, FUN=function(path) {
        verbose && cat(verbose, "Path: ", path);
        unlink(path, recursive=TRUE);
      })
      verbose && exit(verbose);
    }
  }

  # 4. Apply patches?
  if (!is.null(path)) {
    verbose && cat(verbose, "Using patch directory: ", path);
    # (Returns the number of files sourced)
    res <- patchCode(path, verbose=TRUE);
    verbose && exit(verbose);
  } else {
    res <- NULL;
  }
    
  invisible(res);
}) # patchPackage()


############################################################################
# HISTORY:
# 2007-05-10
# o Now patchPackage(..., deleteOld=TRUE) removes not only old patches, but
#   also patches that are older than the package itself.
# 2007-02-28
# o Created.
############################################################################
