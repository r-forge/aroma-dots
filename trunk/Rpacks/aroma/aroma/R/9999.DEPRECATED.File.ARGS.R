#########################################################################/**
# @set "class=File"
# @RdocMethod validateArgumentFileAndPath
#
# @title "Validates the arguments 'file' and 'path'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{A @character string or a @see "R.io::File" object specifying 
#     the file.}
#   \item{path}{A @character string or a @see "R.io::File" object specifying 
#     the path to the file (if the file does not specify that). 
#     If \code{file} is a @see "R.io::File", this argument is ignored.}
#   \item{toWrite}{If @TRUE, specifies that we want to write to the file.}
#   \item{overwrite}{If @TRUE and \code{overwrite} is @TRUE, existing
#     files are overwritten. Otherwise, and Exception is thrown.}
#   \item{mkdir}{If @TRUE, \code{overwrite} is @TRUE, and the path to
#     the file does not exist, it is (recursively) created.}
# }
#
# \value{
#  Returns a @character string of the absolute pathname of the file.
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentFileAndPath", "File", function(this, file=NULL, path=NULL, toWrite=FALSE, overwrite=FALSE, mkdir=toWrite) {
  warning("File$validateArgumentFileAndPath() is deprecated. Please use the Arguments class in R.utils instead.");

  require(R.io);

  if (length(file) > 1) {
    arg <- file;
    argStr <- paste(arg, collapse=", "); 
    throw("The length of argument 'file' is two or more: ", argStr);
  }

  if (length(path) > 1) {
    arg <- path;
    argStr <- paste(arg, collapse=", "); 
    throw("The length of argument 'path' is two or more: ", argStr);
  }

  # Validate 'path'
  if (inherits(path, "File")) {
  } else if (is.character(path)) {
    if (path == "" || path == ".") {
      path <- NULL;   # Current path should be NULL below.
    } else {
      path <- File(path);
    }
  } else if (is.null(path)) {
    # NULL is ok
  } else {
    arg <- path;
    argStr <- paste(arg, collapse=", "); 
    throw("Unknown value of argument 'path': ", argStr);
  }

  # Validate 'file'
  if (inherits(file, "File")) {
    # If file is a File object, the 'path' argument is ignored.
    if (!is.null(path)) {
      path <- NULL;
      warning("Argument 'path' was ignored since argument 'file' is a File object.");
    }
  } else if (is.character(file)) {
    if (is.null(path))
      file <- File(file)
    else
      file <- File(path, file);
  } else if (inherits(file, "connection")) {
    throw("Argument 'file' can not be a connection in this context: ", argStr);
  } else if (is.null(file)) {
    throw("Argument 'file' is NULL.");
  } else {
    arg <- file;
    argStr <- paste(arg, collapse=", "); 
    throw("Unknown value of argument 'file': ", argStr);
  }

  if (toWrite) {
    path <- getParentFile(file);
    if (isExisting(file)) {
      if (!overwrite)
        throw("File already exists and will not be overwritten: ", as.character(file));
    } else if (!isDirectory(path)) {
      # Assert that the path exists.
      if (!mkdir)
        throw("Filepath does not exist: ", as.character(path));
      if (!mkdir(path))
  	throw("Could not create filepath: ", as.character(path));
    }
  }

  file;
}, static=TRUE, protected=TRUE, deprecated=TRUE)


############################################################################
# HISTORY:
# 2005-07-21
# o Made method deprecated.
# 2003-11-18
# o BUG FIX: path=="" and path=="." should be interpreted as path==NULL.
# 2003-10-30
# o Created. See also MicroarrayData.ARGS.R.
############################################################################

