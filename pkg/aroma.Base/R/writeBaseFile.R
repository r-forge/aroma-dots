###########################################################################/**
# @RdocDefault writeBaseFile
#
# @title "Low-level function to write a BASE file structure to a connection or a file"
#
# \description{
#  @get "title".  Note that it is required that the structure is correct;
#  Only minimal validation of structure is done, e.g. the existance of
#  a mandatory header named 'section' is asserted.
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character string filename.}
#   \item{base}{The BASE file @list structure to be written.}
#   \item{verbose}{Either a @logical, a @numeric, or a @see "R.utils::Verbose"
#     object specifying how much verbose/debug information is written to
#     standard output. If a Verbose object, how detailed the information is
#     is specified by the threshold level of the object. If a numeric, the
#     value is used to set the threshold of a new Verbose object. If @TRUE, 
#     the threshold is set to -1 (minimal). If @FALSE, no output is written
#     (and neither is the \link[R.utils:R.utils-package]{R.utils} package required).
#   }
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE, if structure was succesfully written, otherwise @FALSE.
# }
#
# \seealso{
#   @see "readBaseFile".
#   @see "writeBaseFileSection".
#   See \code{stdout()} in \link[base]{showConnections} to write to standard
#   output.
# }
#
# @author
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("writeBaseFile", "default", function(con, base, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (inherits(con, "connection")) {
  } else {
    filename <- as.character(con);

    # Open a file handler
    con <- file(filename, open="w");
    on.exit(close(con), add=TRUE);
  }

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    require(R.utils) || throw("Package not available: R.utils");
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose) {
      require(R.utils) || throw("Package not available: R.utils");
      verbose <- Verbose(threshold=-1);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing BASE file");
  on.exit(verbose && exit(verbose), add=TRUE);

  cat(file=con, "BASEfile\n");

  for (section in base) {
    writeBaseFileSection(con, section, ..., verbose=verbose);
  }

  invisible(TRUE);
}) # writeBaseFile()


############################################################################
# HISTORY: 
# 2005-06-17
# o Tested on several BASE files; seems to work.  Deals also with splitted
#   BASE file structures, i.e. those where 'spots' data tables have been
#   extracted to file by readBaseFile().
# o Created.
############################################################################
