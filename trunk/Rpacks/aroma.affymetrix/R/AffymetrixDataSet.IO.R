###########################################################################/**
# @set "class=AffymetrixDataSet"
# @RdocMethod fromFiles
#
# @title "Defines an AffymetrixDataSet object by searching for Affymetrix files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The path where to search for Affymetrix files.}
#  \item{pattern}{The filename pattern for match files. 
#     If @NULL, filename extensions corresponding to known subclasses
#     of the abstract @see "AffymetrixDataFile" class are search for.}
#  \item{recursive}{If @TRUE, subdirectories are search recursively,
#     otherwise not.}
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "AffymetrixDataSet" object.
# }
#
# \section{Detection of read maps}{
#   If the first data file has a read map specified, then the corresponding
#   map file is searched for and read and set to be the read map for all
#   data files in this data set.  Thus, all data files must share the same
#   map, if any.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "AffymetrixDataSet", function(static, path="cel/", pattern=NULL, recursive=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Get the Class object of the AffymetrixDataFile class corresponding to
  # the static instance 'static'.
  className <- gsub("set$", "File", class(static)[1]);
  clazz <- Class$forName(className);

  # Argument 'pattern':
  if (is.null(pattern)) {
    # Build regular expression pattern from known subclasses
    classNames <- getKnownSubclasses(clazz);

    # File name extensions to search for
    pattern <- gsub("Data", "([A-Z][a-zA-Z0-9]*)", className);
    pattern <- paste("^", pattern, "$", sep="");

    exts <- c();
    for (name in classNames) {
      ext <- gsub(pattern, "\\1", name);
      ext <- tolower(ext);
      exts <- c(exts, tolower(ext), toupper(ext));
    }
    exts <- unique(exts);
    pattern <- paste("[.](", paste(exts, collapse="|"), ")$", sep="");
  } else {
    pattern <- Arguments$getRegularExpression(pattern);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create AffymetrixDataFile object from the matching files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for Affymetrix data files files
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, 
                                                  recursive=recursive, ...);
  if (length(pathnames) == 0)
    throw("No data files found: ", path);

  dataFiles <- list();
  dfStatic <- getStaticInstance(clazz);
  for (kk in seq(along=pathnames)) {
    df <- fromFile(dfStatic, pathnames[kk]);
    dataFiles[[kk]] <- df;
    if (kk == 1) {
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }

  # Create the data-set object
  dataset <- newInstance(static, dataFiles);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify map type
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(dataFiles) > 1) {
    mapType <- getMapType(dataFiles[[1]]);
    if (!is.null(mapType)) {
      apdMap <- ApdMap$fromMapType(mapType);
      setApdMap(dataset, apdMap);
    }
  }

  dataset;
}, static=TRUE)


############################################################################
# HISTORY:
# 2006-05-30
# o Made fromFiles() generic to any subclass.
############################################################################

