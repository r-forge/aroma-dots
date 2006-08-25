###########################################################################/**
# @set "class=AffymetrixDataFile"
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixDataFile object from a file"
#
# \description{
#  @get "title" by calling the same static method of all known subclasses
#  "bottom up" until the data file is "accepted".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of to the file.}
#   \item{path}{The path to the file.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.checkArgs}{Internal.}
# }
#
# \value{
#  Returns an instance of a subclass of @see "AffymetrixDataFile".  
#  If the data file is not recognized by any of the subclasses, an
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
setMethodS3("fromFile", "AffymetrixDataFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
  } else {
    pathname <- filename;
  }

  # Get all known subclasses (bottom up)
  clazz <- Class$forName(class(static)[1]);
  knownSubclasses <- rev(getKnownSubclasses(clazz));
  for (className in knownSubclasses) {
    verbose && enter(verbose, "Trying to define object from file using class ", className);
    clazz <- Class$forName(className);

    # Try reading the file using the static fromFile() method of each class
    static <- getStaticInstance(clazz);
    tryCatch({
      res <- fromFile(static, filename=pathname, .checkArgs=FALSE);
      return(res);
    }, error = function(ex) {
      verbose && cat(verbose, "Failed. Reason: ", ex$message);
    })
    verbose && exit(verbose);
  }

  throw("Tried to read SNP data file using all known subclasses of ", 
         class(static)[1], ": ", paste(knownSubclasses, collapse=", "));
}, static=TRUE)


############################################################################
# HISTORY:
# 2006-07-05
# o Added argument 'verbose'.
# 2006-05-30
# o Created generic static fromFile() which tries to call ditto of all
#   subclasses bottom up.  This is similar to the reading approach taken
#   in the aroma package.
############################################################################
