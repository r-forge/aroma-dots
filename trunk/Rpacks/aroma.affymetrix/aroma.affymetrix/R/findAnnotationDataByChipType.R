###########################################################################/**
# @RdocDefault findAnnotationDataByChipType
#
# @title "Locates an annotation data file by its chip type"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{chipType}{A @list of @see "AffymetrixCelSet".}
#   \item{pattern}{A filename pattern to search for.}
#   \item{...}{Arguments passed to @see "affxparser::findFiles".}
#   \item{paths}{A @character @vector of paths to search.
#     If @NULL, default paths are used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("findAnnotationDataByChipType", "default", function(chipType, pattern=chipType, ..., paths=NULL, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Searching for annotation data file");
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Filename pattern: ", pattern);
  verbose && cat(verbose, "Paths (from argument): ", paste(paths, collapse=", "));

  settings <- getOption("aroma.affymetrix settings");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get search paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paths0 <- paths;

  # Get paths to search
  if (is.null(paths)) {
    paths <- settings$annotationData$paths;
    verbose && cat(verbose, "Paths (from settings): ", paste(paths, collapse=", "));
  }

  if (is.null(paths)) {
    paths <- "annotationData";
    verbose && cat(verbose, "Paths (default): ", paste(paths, collapse=", "));
  } else {
    # Split path strings by semicolons.
    paths <- unlist(strsplit(paths, split=";"));
    verbose && cat(verbose, "Paths (updated): ", paste(paths, collapse=", "));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Expand any file system links
  mainChipType <- gsub("[,].*", "", chipType);
  verbose && cat(verbose, "Main chip type: ", mainChipType);
  paths <- file.path(paths, "chipTypes", mainChipType);
  paths <- sapply(paths, FUN=filePath, expandLinks="any");

  verbose && cat(verbose, "Paths (final): ", paste(paths, collapse=", "));

  # Search recursively for all CDF files
  pathname <- doCall("findFiles", pattern=pattern, paths=paths, recursive=TRUE, ...);

##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   # If not found, look for aliased chip types
##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   if (is.null(pathname)) {
##     verbose && cat(verbose, "No matching annotation file found");
## 
##     settings <- getOption("aroma.affymetrix.settings");
##     aliases <- settings$annotationData$aliases;
##     verbose && cat(verbose, "settings$annotationData$aliases:");
##     verbose && print(verbose, aliases);
## 
##     if (!is.null(aliases)) {
##       verbose && enter(verbose, "Looking for aliased chip types");
## 
##       aliases <- as.character(aliases);
##       keys <- names(aliases);
##       if (is.null(keys)) {
##         # Assume format "foo=val1;bar=val2".
##         # Turn into a named character vector.
##         parts <- unlist(strsplit(aliases, split=";"));
##         parts <- strsplit(parts, split="=");
##         keys <- unlist(lapply(parts, FUN=function(x) x[1]), use.names=FALSE);
##         aliases <- unlist(lapply(parts, FUN=function(x) x[2]), use.names=FALSE);
##         names(aliases) <- keys;
##       }
##       verbose && cat(verbose, "Parsed aliases:");
##       verbose && print(verbose, aliases);
## 
##       # Get the aliased chip type
##       realChipType <- aliases[chipType];
##       if (!is.null(realChipType) && !is.na(realChipType)) {
##         verbose && cat(verbose, "Aliased chip type:");
##         verbose && print(verbose, realChipType);
##         
##         # AD HOC: Replace all occurances of <chipType> in the pattern 
##         # with <realChipType>.
##         pattern <- gsub(chipType, realChipType, pattern, fixed=TRUE);
##         verbose && cat(verbose, "Filename pattern (updated): ", pattern);
##         verbose && enter(verbose, "Searching for aliased chip type instead");
##         pathname <- findAnnotationDataByChipType(realChipType, 
##                                      pattern=pattern, ..., paths=paths0);
##         verbose && exit(verbose);
##       } else {
##         verbose && cat(verbose, "No alias available: ", chipType);
##       }
## 
##       verbose && exit(verbose);
##     } # if (!is.null(aliases))
##   }

  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && exit(verbose);

  pathname;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2007-02-21
# o Added Rdoc comments.
# o Added verbose.
# o Added support for aliases.
# o Changed settings$paths$annotationData to settings$annotationData$paths.
# 2007-02-06
# o Created.
############################################################################
