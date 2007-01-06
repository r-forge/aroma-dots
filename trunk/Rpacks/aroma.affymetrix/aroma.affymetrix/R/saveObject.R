###########################################################################/**
# @RdocDefault saveObject
#
# @title "Saves an object to a file or a connection"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{object}{The object to be saved.}
#  \item{file}{A filename or @connection where the object should be saved.
#    If @NULL, the filename will be the hash code of the object plus ".xdr".}
#  \item{path}{Optional path, if \code{file} is a filename.}
#  \item{compress}{If @TRUE, the file is compressed to, otherwise not.}
#  \item{...}{Other arguments accepted by \code{save()} in the base package.}
# }
#
# \value{
#  Returns (invisibly) the pathname or the @connection.
# }
#
# @author
#
# \seealso{
#   @see "loadObject" to load an object from file.
#   @see "digest::digest" for how hash codes are calculated from an object.
#   Internally @see "base::save" is used.
# }
#
# @keyword programming
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("saveObject", "default", function(object, file=NULL, path=NULL, compress=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file':
  if (is.null(file)) {
    require("digest") || throw("Package 'digest' not loaded.");
    file <- digest(as.list(object));  # Might be slow.
    file <- sprintf("%s.xdr", file);
  } 
  if (!inherits(file, "connection")) {
    file <- filePath(path, file, expandLinks="any");
  }

  saveLoadReference <- object;
  base::save(saveLoadReference, file=file, ..., compress=compress, ascii=FALSE);

  invisible(file);
}) # saveObject()





###########################################################################/**
# @RdocDefault loadObject
#
# @title "Method to load object from a file or a connection"
#
# \description{
#   @get "title", which previously have been saved using @see "saveObject".
# }
#
# @synopsis
#
# \arguments{
#  \item{file}{A filename or @connection to read the object from.}
#  \item{path}{The path where the file exists.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the save object.
# }
#
# \details{
#   The main difference from this method and @see "base::load" in the
#   \pkg{base} package, is that this one returns the object read rather
#   than storing it in the global environment by its default name.
#   This makes it possible to load objects back using any variable name.
# }
#
# @author
#
# \seealso{
#   @see "saveObject" to save an object to file.
#   Internally @see "base::load" is used.
# }
#
# @keyword programming
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("loadObject", "default", function(file, path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file':
  if (!inherits(file, "connection")) {
    file <- filePath(path, file, expandLinks="any");
  }
 
  # load.default() recognized gzip'ed files too.
  vars <- base::load(file=file);

  if (!"saveLoadReference" %in% vars)
    throw("The file was not save by saveObject(): ", file);

  saveLoadReference;
}) # loadObject()



##############################################################################
# HISTORY:
# 2006-11-24
# o Created from Object.R in the R.oo package. This will probably be moved
#   to either R.oo or R.utils later.
##############################################################################
