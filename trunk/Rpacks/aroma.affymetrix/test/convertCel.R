#########################################################################/**
# @RdocFunction convertCel
#
# @title "Converts a CEL into the same CEL but with another format"
#
# @synopsis 
# 
# \description{
#   @get "title".  
#   Currently only CEL files in version 4 (binary/XDA) can be written.
#   However, any input format is recognized.
# }
# 
# \arguments{
#   \item{filename}{The pathname of the original CEL file.}
#   \item{outFilename}{The pathname of the destination CEL file.
#     If the same as the source file, an exception is thrown.}
#   \item{version}{The version of the output file format.}
#   \item{...}{Not used.}
#   \item{verbose}{If @TRUE, extra details are written while processing.}
# }
# 
# \value{
#   Returns (invisibly) @TRUE if a new CEL was generated, otherwise @FALSE.
# }
#
# \section{Benchmarking of ASCII and binary CELs}{
#   Binary CELs are much faster to read than ASCII CELs.  Here are some
#   example for reading complete CELs (the differnce is even larger when
#   reading CELs in subsets):
#   \itemize{
#     \item To do
#   }
# }
#
# @examples "../incl/convertCel.Rex"
#
# @author
#
# \seealso{
#   @see "createCel".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
convertCel <- function(filename, outFilename, version="4", force=FALSE, ..., .validate=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  # Expand any '~' in the pathname.
  filename <- file.path(dirname(filename), basename(filename));
  if (!file.exists(filename)) {
    stop("Cannot open CEL file. File does not exist: ", filename);
  }

  # Argument 'outFilename':
  # Expand any '~' in the pathname.
  outFilename <- file.path(dirname(outFilename), basename(outFilename));
  if (identical(outFilename, filename)) {
    stop("Cannot convert CEL file. Destination is identical the the source pathname: ", filename);
  }

  # Argument 'version':
  version <- as.character(version);
  if (version == "4") {
  } else {
    stop("Cannot convert CEL. Currently only version 4 (binary/XDA) can be written: ", version);
  }

  # Argument 'verbose':
  verbose <- as.integer(verbose);

  verbose2 <- verbose-1;
  if (verbose2 < 0) verbose2 <- 0;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read source CEL
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (verbose)
    cat("Reading CEL file...\n");
  cel <- readCel(filename, readHeader=TRUE, readXY=TRUE, readIntensities=TRUE, readStdvs=TRUE, readPixels=TRUE, readOutliers=FALSE, readMasked=FALSE);
  if (verbose)
    cat("Reading CEL file...done\n");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Creating new CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (verbose)
    cat("Creating empty CEL file...\n");
  cel$header$version <- 4;
  pathname <- createCel(outFilename, header=cel$header, overwrite=FALSE, verbose=verbose2);
  if (verbose)
    cat("Creating empty CEL file...done\n");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update destination CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (verbose)
    cat("Updating CEL file...\n");
  indices <- cel$X + cel$Y*cel$header$cols;
  updateCel(outFilename, indices=indices, intensities=cel, verbose=verbose2);
  if (verbose)
    cat("Updating CEL file...done\n");

  if (.validate) {
    compareCels(filename, outFilename, verbose=verbose2);
  }

  invisible(pathname);
} # convertCel()


############################################################################
# HISTORY:
# 2007-01-03
# o Created.
############################################################################
