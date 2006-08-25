###########################################################################/**
# @RdocClass AffymetrixApdFile
#
# @title "The AffymetrixApdFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixApdFile object represents a single APD file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AffymetrixApdFile", function(...) {
  extend(AffymetrixDataFile(...), "AffymetrixApdFile")
})

setMethodS3("getChipType", "AffymetrixApdFile", function(this, ...) {
  pathname <- getPathname(this);
  readApdHeader(pathname)$chipType;
})

setMethodS3("getFileType", "AffymetrixApdFile", function(this, ...) {
  "apd";
})

setMethodS3("nbrOfProbes", "AffymetrixApdFile", function(this, ...) {
  pathname <- getPathname(this);
  # Affymetrix probe data is stored in a file vector.
  apd <- FileVector(pathname);
  on.exit(close(apd));
  length(apd);
})

setMethodS3("getMapType", "AffymetrixApdFile", function(this, ...) {
    pathname <- getPathname(this);
    header <- readApdHeader(pathname);
    header$mapType;
}, protected=TRUE)


setMethodS3("readIntensities", "AffymetrixApdFile", function(this, probes=NULL, ...) {
  pathname <- getPathname(this);
  readApd(pathname, indices=probes, ..., name="intensities")$intensities;
})

setMethodS3("readUnits", "AffymetrixApdFile", function(this, ...) {
  pathname <- getPathname(this);
  readApdUnits(pathname, ...);
})

############################################################################
# HISTORY:
# 2006-07-21
# o Added readUnits().
# 2006-03-18
# o Made probe indices one-based.
# 2006-03-02
# o Created.
############################################################################
