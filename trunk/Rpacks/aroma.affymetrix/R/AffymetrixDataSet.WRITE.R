###########################################################################/**
# @set "class=AffymetrixDataSet"
# @RdocMethod writeApd
#
# @title "Generates APD files for each of the data files in this data set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mapType}{The type of read map for the generated APD file.
#     If @NULL, no remapping of the cell indices is done.
#     If \code{"asChipType"}, the map type is the same as the chip type
#     of the CEL file.  
#     If any other @character string, it sets the map type to that string.
#     Note that there must be a APD map file with this type that can
#     be found by @see "aroma.apd::findApdMap".
#   }
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a new \code{AffymetrixDataSet} instance.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("writeApd", "AffymetrixDataSet", function(this, apdMap="byChipType", ...) {
  # Argument 'apdMap':
  if (is.null(apdMap)) {
  } else if (identical(apdMap, "byChipType")) {
    chipType <- getChipType(this);

    # Create an optimal read map
    apdMap <- ApdMap$fromCdf(chipType=chipType);

    # Write map to current directory if non existing
    if (is.null(findApdMap(apdMap))) {
      write(apdMap);
    }
  } else if (!inherits(apdMap, "ApdMap")) {
    throw("Argument 'apdMap' is not an ApdMap object: ", class(apdMap)[1]);
  }

  # For each of the data files in this set, call writeApd().
  dataFiles <- lapply(this, FUN=writeApd, apdMap=apdMap, ...);

  # Creates a new AffymetrixDataSet object for the APD data files
  res <- newInstance(this, dataFiles);

  res;
})



###########################################################################/**
# @RdocMethod write
#
# @title "Writes the data files in this data set to file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{format}{A @character string specifying the file format to be written.}
#  \item{...}{Arguments passed to the internal write function.}
# }
#
# \value{
#   Returns a new \code{AffymetrixDataSet} instance.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("write", "AffymetrixDataSet", function(this, format="apd", ...) {
  # Argument 'format':
  format <- match.arg(format);

  # Find the write function
  name <- paste("write", capitalize(format), sep="");
  fcn <- get(name, mode="function");

  # For each of the data files in this set, call the write function.
  dataFiles <- lapply(this, FUN=fcn, ...);

  fcn(this, ...);
})



############################################################################
# HISTORY:
# 2006-05-15
# o Extracted from AffymetrixDataSet.R.
# 2006-03-04
# o Added writeApd().
############################################################################
