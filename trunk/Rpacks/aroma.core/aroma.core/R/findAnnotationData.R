###########################################################################/**
# @RdocDefault findAnnotationData
#
# @title "Locates an annotation data file"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{name}{A @character string.}
#   \item{tags}{Optional @character string.}
#   \item{pattern}{A filename pattern to search for.}
#   \item{private}{If @FALSE, files and directories starting 
#     with a periods are ignored.}
#   \item{...}{Arguments passed to @see "affxparser::findFiles".}
#   \item{paths}{A @character @vector of paths to search.
#     If @NULL, default paths are used.}
#   \item{set}{A @character string specifying what type of annotation 
#     to search for.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("findAnnotationData", "default", function(name, tags=NULL, set, pattern=NULL, private=FALSE, ..., paths=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'set':
  set <- Arguments$getCharacter(set);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fullname <- paste(c(name, tags), collapse=",");
  if (is.null(pattern)) {
    pattern <- fullname;
  }

  name <- gsub("[,].*", "", fullname);

  verbose && enter(verbose, "Searching for annotation data file(s)");
  verbose && cat(verbose, "Name: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
  verbose && cat(verbose, "Set: ", set);
  verbose && cat(verbose, "Filename pattern: ", pattern);
  verbose && cat(verbose, "Paths (from argument): ", paste(paths, collapse=", "));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get search paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paths0 <- paths;

  if (is.null(paths)) {
    paths <- "annotationData";
    verbose && cat(verbose, "Paths (default): ", paste(paths, collapse=", "));
  } else {
    # Split path strings by semicolons.
    paths <- unlist(strsplit(paths, split=";"));
    verbose && cat(verbose, "Paths (updated): ", paste(paths, collapse=", "));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/<set>/<genome>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Expand any file system links
  paths <- file.path(paths, set, name);
  paths <- sapply(paths, FUN=filePath, expandLinks="any");

  verbose && cat(verbose, "Paths (final): ", paste(paths, collapse=", "));

  # Search recursively for all files
  args <- list(...);
  args$pattern <- pattern;
  args$allFiles <- private;
  args$paths <- paths;
  args$recursive <- TRUE;
  verbose && cat(verbose, "Arguments to findFiles:");
  verbose && str(verbose, args);
  pathname <- do.call("findFiles", args=args);

  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && exit(verbose);

  pathname;
}, protected=TRUE)  # findAnnotationData()

############################################################################
# HISTORY:
# 2008-05-09
# o Removed the option to specify the annotation data path by the option
#   'aroma.affymetrix.settings'.
# 2008-02-14
# o Added argument 'private=TRUE' to findAnnotationData().
# 2007-09-15
# o Created from findAnnotationDataByChipType.R.
############################################################################
