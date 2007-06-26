###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod getSiblings
#
# @title "Gets the all data sets that refers to the same samples a this one"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{notSelf}{If @TRUE, the current data set is not returned, just
#      its siblings.}
#   \item{...}{Additional arguments passed to @seemethods "fromFiles".}
# }
#
# \value{
#  Returns a named @list of @see "AffymetrixCelSet" objects.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getSiblings", "AffymetrixCelSet", function(this, notSelf=FALSE, ...) {
  # Scan parent directory for all possible data sets.
  # /path/to/<data set name>/<chip type>/
  path <- getPath(this);

  # /path/to/<data set name>/
  parent <- getParent(path);

  # Scan for all directories
  paths <- list.files(path=parent, full.names=TRUE);
  isDir <- sapply(paths, isDirectory);
  paths <- paths[isDir];

  if (notSelf) {
    paths <- paths[paths != getPath(this)];
  }

  # Keep only directories containing CEL files
  hasCELs <- sapply(paths, FUN=function(path) {
    files <- list.files(path=path, pattern="[.](cel|CEL)$");
    (length(files) > 0);
  })
  paths <- paths[hasCELs];

  sets <- vector("list", length(paths));
  names(sets) <- basename(paths);
  for (kk in seq(along=paths)) {
    path <- paths[kk];
    if (!notSelf && identical(path, getPath(this))) {
      sets[[kk]] <- this;
    } else {
      sets[[kk]] <- fromFiles(this, path=path, ...);
    }
  }
  
  sets;
}, private=TRUE)




############################################################################
# HISTORY:
# 2007-06-25
# o Now '...' are passed to fromFiles() in getSiblings().
# 2007-06-12
# o BUG FIX: getSiblings() for AffymetrixCelSet was broken.
# o Moved into its own *.R file.
# 2007-04-06
# o BUG FIX: fromFiles() of AffymetrixCelSet would give error "Exception: 
#   Pathname not found: annotationData/samples" if that directory was 
#   missing.  Now it is instead created.
# 2007-04-03
# o Now fromFiles() verifies that the set CDF is compatible with the CEL
#   files, otherwise the CDF of the first CEL file is used.
# 2007-04-02
# o Now sample annotation files are searched for in annotationData/samples/.
# 2007-03-28
# o Added argument 'cache=TRUE' to getUnitIntensities() and readUnits().
# 2007-03-24
# o BUG FIX: clearCache() did not clear the .readUnitsCache field.
# 2007-03-16
# o BUG FIX: getAverageFile() of AffymetrixCelSet would average the wrong
#   set of cells if argument 'indices' was different from NULL.
# 2007-03-06
# o Now attributes are set from SAF files in fromFiles().
# 2007-02-22
# o Fixed the warning about "'tzone' attributes are inconsistent". See
#   code of as.character() for explanation.
# o Now fromFiles() accepts argument 'chipType' to override any chip type
#   specified in the CEL headers. This is useful in case different CEL files
#   refers to different chip types, which can be the case for mixed 
#   generations of CEL files.  Also added a scan of chip types.
# 2007-02-14
# o Added test for correct directory structure to fromFiles().  This will
#   enforce users to use the correct structure so that for instance the
#   name of the data set is correctly inferred.
# o Added argument 'onDuplicates' to fromFiles() so it is possible to 
#   exclude duplicated CEL files.
# 2007-02-06
# o Added static method findByName() and fromName().
# 2007-01-15
# o Added 'classes=class(this)' to all "digest" keys.
# 2007-01-07
# o BUG FIX: In KS's update of getAverageFile() to support averaging
#   over other fields than intensities, argument 'indices' was missing
#   in the readCel() call making the function fail when processed chunk
#   by chunk. /HB
# 2007-01-05
# o Removed getSampleNames().
# 2006-12-01
# o Now as.character() reports the range of CEL header timestamps.
# 2006-11-07
# o Now getAverageFile() uses rowMedians() of R.native if available, 
#   otherwise a local version utilizing apply(). Same for rowMads().
# 2006-10-24
# o Added getAverageLog() and getAverageAsinh().
# o Added transforms and anti-transforms g() and h() to getAverageFile().
# o Changed the defaults from mean to median, and sd to mad for 
#   getAverageFile().
# o Added Rdoc comments to getAverageFile().
# 2006-10-10
# o Renamed rma and gcrma to rmaSummary and gcrmaSummary, to avoid clash
#   with existing functions.
# o Added gcrma() wrapper function.
# o Added rma() wrapper function.
# o Fixed a bug in getData() - default for argument "fields" contained "xy",
#   which is not a valid field (x, y are separate).
# 2006-10-02
# o Added getData().  Now getIntensities() works again and is just a wrapper
#   to getData().
# 2006-09-18
# o Now references to all requested average files are cached so it can
#   return the same object instead of creating a new one each time.
# 2006-09-16
# o Added getSiblings() to easily get other data sets for the same
#   samples.
# 2006-09-14
# o Added a read-buffer cache to readUnits() and getUnitIntensities().
# 2006-08-27
# o Added getAverageFile().
# 2006-08-26
# o Now getName() of a CEL set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Now AffymetrixCelSet inherits from AffymetrixFileSet.
# 2006-08-11
# o Added clearCache() which also clears the cache of all data file object.
# 2006-05-16
# o Redefined "[" to extract arrays.
# 2006-04-13
# o Added Rdoc comments for all methods.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-30
# o Updated to new aroma.apd.
# 2006-03-18
# o Added argument 'subset' to calcAvgCellSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfCells() returns the number of cells for the first file only.
# o Now the fromFiles(static, ...) creates an object of the same class as 
#   the static object.
# 2006-03-04
# o Added mapping functions.
# o Added writeApd().
# 2006-03-03
# o Added lapply().
# 2006-03-02
# o Updated to deal with AffymetrixDataFile object instead of CEL directly.
# 2006-02-21
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
