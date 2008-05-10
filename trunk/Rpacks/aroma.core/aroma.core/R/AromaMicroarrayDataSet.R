###########################################################################/**
# @RdocClass AromaMicroarrayDataSet
#
# @title "The AromaMicroarrayDataSet class"
#
# \description{
#  @classhierarchy
#
#  An AromaMicroarrayDataSet object represents a set of 
#  @see "AromaMicroarrayDataFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AromaMicroarrayDataFile":s.}
#   \item{...}{Arguments passed to @see "GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaMicroarrayDataSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    reqFileClass <- "AromaMicroarrayDataFile";
    base::lapply(files, FUN=function(df) {
      if (!inherits(df, reqFileClass))
        throw("Argument 'files' contains a non-", reqFileClass, 
                                                      " object: ", class(df));
    })
  } else if (inherits(files, "AromaMicroarrayDataSet")) {
    return(as.AromaMicroarrayDataSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(GenericDataFileSet(files=files, ...), "AromaMicroarrayDataSet");
})



setMethodS3("clearCache", "AromaMicroarrayDataSet", function(this, ...) {
  # Clear the cache of all files
  lapply(this, clearCache);

  # Then for this object
  NextMethod("clearCache", this);
})




###########################################################################/**
# @RdocMethod as.AromaMicroarrayDataSet
# @alias as.AromaMicroarrayDataSet.list
# @alias as.AromaMicroarrayDataSet.default
#
# @title "Coerce an object to an AromaMicroarrayDataSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "AromaMicroarrayDataSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AromaMicroarrayDataSet", "AromaMicroarrayDataSet", function(object, ...) {
  object;
})

setMethodS3("as.AromaMicroarrayDataSet", "list", function(object, ...) {
  AromaMicroarrayDataSet(object, ...);
})

setMethodS3("as.AromaMicroarrayDataSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AromaMicroarrayDataSet object: ", 
                                                                mode(object));
})



###########################################################################/**
# @RdocMethod fromFiles
#
# @title "Defines an AromaMicroarrayDataSet object by searching for data files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The path where to search for data files.}
#  \item{...}{Optional arguments passed to the same method in the superclass.}
#  \item{.validate}{If @TRUE, the data set it validated.}
# }
#
# \value{
#   Returns an @see "AromaMicroarrayDataSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "AromaMicroarrayDataSet", function(static, ..., fileClass="AromaMicroarrayDataFile", .validate=TRUE) {
  # NextMethod() does not work here
  res <- fromFiles.GenericDataFileSet(static, ..., fileClass=fileClass);

  if (.validate) {
    chipTypes <- lapply(res, FUN=getChipType);
    chipTypes <- unique(chipTypes);
    if (length(chipTypes) > 1) {
      throw("The located ", class(res)[1], " contains files with different chip types: ", paste(chipTypes, collapse=", "));
    }
  }

  res;
}, static=TRUE)


setMethodS3("getPlatform", "AromaMicroarrayDataSet", function(this, ...) {
  file <- getFile(this, 1);
  getPlatform(file);
})


setMethodS3("getChipType", "AromaMicroarrayDataSet", function(this, ...) {
  file <- getFile(this, 1);
  getChipType(file);
})


############################################################################
# HISTORY:
# 2008-05-09
# o Created.
############################################################################
