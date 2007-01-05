###########################################################################/**
# @RdocDefault findCelSet
#
# @title "Searches for one or several CEL sets"
#
# \description{
#   @get "title" that matches certain criterias.
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A regular expression pattern to match the name of the
#    \emph{data set} subdirectory.}
#  \item{chipType}{A regular expression pattern to match the name of
#    the \emph{chip type} subdirectory.}
#  \item{paths}{A @character @vector of paths to be searched.
#    If @NULL certain default directories are scanned. See below.}
#  \item{minCount}{The minimum number of CEL files in a directory for
#    it to be considered a data set.}
#  \item{...}{Not used.}
#  \item{firstOnly}{If @TRUE, the first possible match is returned, 
#    otherwise all matches are returned.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns @character @vector of full pathname to matching data set.
# }
#
# \details{
#   This function searches for data sets with a pathname pattern
#   \preformatted{
#      path/to/<data set name>/chip_data/<chip type>/
#   }
#   of a directory containing CEL files.
# }
#
# \section{Default search paths}{
#   Note, the current directory is always searched at the beginning.
#   This provides an easy way to override other files in the search path.
#   If \code{paths} is @NULL, then a set of default paths are searched.
#   The default search path is consituted of:
#   \enumerate{
#    \item \code{"."}
#    \item \code{getOption("AFFX_DATA_PATH")}
#    \item \code{Sys.getenv("AFFX_DATA_PATH")}
#    \item \code{"data/"}
#   }
#
#   One of the easiest ways to set system variables for \R is to
#   set them in an \code{.Renviron} file, e.g.
#   \preformatted{
#     # affxparser: Set default CDF path
#     AFFX_DATA_PATH=${AFFX_DATA_PATH};M:/data
#   }
#   See @see "base::Startup" for more details.
# }
#
# @examples "../incl/findCelSet.Rex"
#
# @author
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("findCelSet", "default", function(name=NULL, chipType=NULL, paths=NULL, minCount=1, ..., firstOnly=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # From affxparser::findFiles() 
  splitPaths <- function(paths, ...) {
    if (length(paths) == 0)
      return(NULL);
    # If in format "path1; path2;path3", split it to multiple strings.
    paths <- unlist(strsplit(paths, split=";"));
    paths <- gsub("[ \t]*$", "", gsub("^[ \t]*", "", paths));
    paths <- paths[nchar(paths) > 0];
    if (length(paths) == 0)
      return(NULL);
    paths;
  }

  # If R.utils::filePath() is missing
  hasRUtils <- require(R.utils);
  if (!hasRUtils) {
    filePath <- function(..., expandLinks="any") {
      file.path(...);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paths':
  if (is.null(paths)) {
    paths <- paste(".", "data/", getOption("AFFX_DATA_PATH"), 
           Sys.getenv("AFFX_DATA_PATH"), sep = ";", collapse = ";");
  }
  paths <- splitPaths(paths);
  paths <- c(".", paths);
  paths <- unique(paths);

  # Argument 'name':
  if (!is.null(name)) {
    if (hasRUtils)
      name <- Arguments$getRegularExpression(name);
  }

  # Argument 'chipType':
  if (!is.null(chipType)) {
    if (hasRUtils)
      chipType <- Arguments$getRegularExpression(chipType);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  allPaths <- c();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan each search path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Scanning ", length(paths), " paths");
  verbose && print(verbose, paths);
  for (path in paths) {
    path <- filePath(path, expandLinks="any");
    verbose && print(verbose, path);
    if (!isDirectory(path))
      next;

    # Retrieve all directories at this path
    dirs <- list.files(path=path, pattern=name, full.names=TRUE);
    if (length(dirs) == 0)
      next;
    dirs <- dirs[sapply(dirs, FUN=isDirectory)];
    if (length(dirs) == 0)
      next;
    verbose && cat(verbose, "Directories:");
    verbose && str(verbose, dirs);
  
    # Keep only directories with a chip_data/ subdirectory
    dirs <- dirs[sapply(dirs, FUN=function(dir) {
      subpath <- filePath(dir, "chip_data", expandLinks="any");
      isDirectory(subpath);
    })]
    if (length(dirs) == 0)
      next;

    verbose && cat(verbose, "Directories with a chip_data/ subdirectory:");
    verbose && str(verbose, dirs);

    # The path to all such
    dirs <- sapply(dirs, FUN=filePath, "chip_data", expandLinks="any");
  
    # Select by chip type(s)
    dirs <- sapply(dirs, FUN=function(dir) {
      list.files(path=dir, pattern=chipType, full.names=TRUE);
    })
    dirs <- dirs[sapply(dirs, function(x) (length(x) > 0))]
    dirs <- unname(dirs);
    dirs <- unlist(dirs, use.names=FALSE);
    if (length(dirs) == 0)
      next;

    # Count the number of CEL files
    if (minCount > 0) {
      dirs <- dirs[sapply(dirs, function(path) {
        files <- list.files(pattern="[.](c|C)(e|E)(l|L)$", path=path);
        (length(files) >= minCount);
      })];
    }

    if (firstOnly && length(dirs) > 0)
      return(dirs[1]);  

    allPaths <- c(allPaths, dirs);
  }
  verbose && exit(verbose);

  allPaths;
}, private=TRUE) # findCelSet()


############################################################################
# HISTORY:
# 2006-09-15
# o Created.  This will simplify retrieving Affymetrix data sets, but it 
#   will also standardize how data sets are stored.
############################################################################
