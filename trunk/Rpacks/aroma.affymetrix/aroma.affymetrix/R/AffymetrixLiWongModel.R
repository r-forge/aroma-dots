###########################################################################/**
# @RdocClass AffymetrixLiWongModel
#
# @title "The AffymetrixLiWongModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Li \& Wong (2001) model,
#  see @see "AffymetrixLiWongModel".
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet" object.}
#   \item{path}{The @character string specifying the path to the directory
#      to contain the parameter-estimate files.}
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \references{
#   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
#   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
# }
#*/###########################################################################
setConstructorS3("AffymetrixLiWongModel", function(dataSet=NULL, path=filePath("modelLiWong", getChipType(getCdf(dataSet))), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affy) || throw("Package 'affy' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(dataSet)) {
    # A work-around for the fact that getCdf(NULL) is not work.
    path=NULL;
  }

  extend(ProbeLevelModel(dataSet=dataSet, path=path, ...), "AffymetrixLiWongModel")
})


setMethodS3("getProbeAffinityClass", "AffymetrixLiWongModel", function(static, ...) {
  LiWongProbeAffinityFile;
}, static=TRUE, protected=TRUE)


setMethodS3("getFitFunction", "AffymetrixLiWongModel", function(static, ...) {
  liWong <- function(y, ...) {
    affy::fit.li.wong(t(y));
  }

  liWong;
}, static=TRUE, protected=TRUE)



setMethodS3("getChipEffects", "AffymetrixLiWongModel", function(this, ..., verbose=FALSE) {
  chipFiles <- this$.chipFiles;
  if (!is.null(chipFiles))
    return(chipFiles);

  # For each of the data files, create a file to store the estimates in
  path <- getPath(this);
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  unitSizes <- getUnitSizes(cdf);
  unitOffsets <- cumsum(unitSizes) - unitSizes[1] + 1;
  n <- sum(unitSizes);

  chipFiles <- list();
  for (kk in seq(ds)) {
    df <- ds[[kk]];
    filename <- sprintf("%s-liwong.apd", getName(df));
    filename <- filePath(path, filename);
    if (!isFile(filename)) {
      X <- FileFloatVector(filename, length=n);
      X <- FileFloatVector(length=n, appendTo=X);
      X <- FileByteVector(length=n, appendTo=X);
      close(X);
      rm(X);
    }
    set <- AbstractFileArray$fromFile(filename);
    # We have to close the parameter files, because Windows can only
    # handle ~512 open connections, and we might have more arrays.
    # Instead we have to open and close the connections each time we
    # read data.
    close(set[[1]]);
    chipFiles[[kk]] <- set;
  } # for (kk in ...)

  this$.chipFiles <- chipFiles;

  chipFiles;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added getProbeAffinities() and the corrsponding cached fields.
# o Now fit() does not re-read data just updated.
# 2006-08-19
# o After all the bug fixes in updateCel() I think this function finally
#   works correctly.
# 2006-08-17
# o Created.
############################################################################
