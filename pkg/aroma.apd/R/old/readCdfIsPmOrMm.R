#########################################################################/**
# @RdocFunction readCdfIsPmOrMm
#
# @title "Checks if cells each group of each unit are mis-match probes"
#
# @synopsis 
# 
# \description{
# }
# 
# \arguments{
#   \item{filename}{The filename of the CDF file.}
#   \item{units}{An @integer @vector of (zero-based) unit indices 
#     specifying which units to be read.  If @NULL, all units are read.}
#   \item{verbose}{An @integer specifying the verbose level. If 0, the 
#     file is parsed quietly.  The higher numbers, the more details.}
# }
# 
# \value{
#   A named @list of named @logical vectors.  The name of the list elements
#   are unit names and the names of the logical vector are group names.
# }
# 
# @examples "../incl/readCdfIsPmOrMm.Rex"
#
# @author
# 
# \seealso{
# }
# 
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
readCdfIsPmOrMm <- function(filename, units=NULL, verbose=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename':
  filename <- file.path(dirname(filename), basename(filename));
  if (!file.exists(filename))
    stop("File not found: ", filename);

  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- as.integer(units);
    if (any(units < 0))
      stop("Argument 'units' contains negative indices.");
  } else {
    stop("Argument 'units' must be numeric or NULL: ", class(units)[1]);
  }

  # Argument 'units':
  if (length(verbose) != 1)
    stop("Argument 'units' must be a single integer.");
  verbose <- as.integer(verbose);
  if (!is.finite(verbose))
    stop("Argument 'units' must be an integer: ", verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read the CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  .Call("R_affx_cdf_isPmOrMm", filename, units, verbose, PACKAGE="affxparser");
}


############################################################################
# HISTORY:
# 2006-01-11
# o Created. /HB
############################################################################
