###########################################################################/**
# @set "class=AffymetrixApdFile"
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixApdFile object from a file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixApdFile" object.  
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AffymetrixApdFile", function(static, filename, path=NULL, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  # Try to read the header
  header <- readApdHeader(pathname);

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE)




############################################################################
# HISTORY:
# 2006-05-30
# o Added fromFile().
############################################################################
