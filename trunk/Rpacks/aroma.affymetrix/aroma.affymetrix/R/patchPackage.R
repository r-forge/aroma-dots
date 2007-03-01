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

    verbose && enter(verbose, "Root path: ", path);
  
    # Search for patch directories
    pattern <- "20[0-9][0-9][01][0-9][0-3][0-9](|[a-z])";
    paths <- list.files(path=rootPath, pattern=pattern, full.names=TRUE);
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
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Argument 'deleteOld':
  deleteOld <- Arguments$getVerbose(deleteOld);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Patching package ", pkgName);

  # 1. Scan for patch directories
  paths <- lapply(paths, FUN=function(path) {
    findPatchDirectories(pkgName, path=path, verbose=verbose);
  })
  paths <- unlist(paths, use.names=FALSE);

  # Nothing do to?
  if (length(paths) == 0) {
    return(invisible(0));
  }

  # 2. Keep only the most recent one
  o <- order(basename(paths), decreasing=TRUE);
  paths <- paths[o];
  path <- paths[1];

  # 3. Remove all older patches
  paths <- paths[-1];
  if (deleteOld && length(paths) > 0) {
    verbose && enter(verbose, "Deleting deprecated patch directories");
    res <- sapply(paths, FUN=function(path) {
      verbose && cat(verbose, "Path: ", path);
      unlink(path, recursive=TRUE);
    })
    verbose && exit(verbose);
  }

  # 4. Apply patches  
  verbose && cat(verbose, "Using patch directory: ", path);
  # (Returns the number of files sourced)
  res <- patchCode(path, verbose=TRUE);
  verbose && exit(verbose);
    
  invisible(res);
}) # patchPackage()


############################################################################
# HISTORY:
# 2007-02-28
# o Created.
############################################################################
