###########################################################################/**
# @RdocDefault writeBaseFileSection
#
# @title "Low-level function to write a BASE file section to a connection or a file"
#
# \description{
#  @get "title". This a supportive function to writeBaseFile().
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character string filename.}
#   \item{section}{The BASE file section @list structure to be written.}
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
#   @see "writeBaseFile".
# }
#
# @author
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("writeBaseFileSection", "default", function(con, section, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function definitions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeBaseFileHeaders <- function(con, section, ...) {
    headers <- section$headers;

    if (length(headers) == 0)
      throw("No headers found in BASE section.");

    if (!"section" %in% names(headers)) {
      throw("Mandatory header 'section' is missing: ", 
                                    paste(names(headers), collapse=", "));
    }

    verbose && enter(verbose, "Writing headers");
    on.exit(verbose && exit(verbose), add=TRUE);

    # Remove any 'dataFiles' headers.
    headers[["dataFiles"]] <- NULL;

    # Remove any 'count' headers from 'spots' section, because they will
    # be written when the data table is written.
    if (headers[["section"]] == "spots")
      headers[["count"]] <- NULL;

    # Write headers
    for (kk in seq(length(headers))) {
      key <- names(headers)[kk];
      values <- headers[[kk]];
      values <- paste(values, collapse="\t");
      cat(file=con, key, "\t", values, "\n", sep="");
    }

    TRUE;
  } # writeBaseFileHeaders()
  
  
  writeBaseFileData <- function(con, section, ...) {
    headers <- section$headers;
    data <- section$data;

    verbose && enter(verbose, "Writing data");
    on.exit(verbose && exit(verbose), add=TRUE);

    verbose && cat(verbose, "Re-constructing data table.")
    dataFiles <- headers[["dataFiles"]];
    if (!is.null(dataFiles)) {
      # First, re-create data table
      assays <- headers[["assays"]];
      columnFields <- setdiff(headers[["columns"]], "assayData");
      assayFields <- headers[["assayFields"]];
      data <- NULL;
      for (dataFile in dataFiles) {
        # Read all data...
        tmp <- read.table(file=dataFile, header=TRUE, sep="\t");
        # ...but keep only wanted fields
        if (is.null(data)) {
          data <- tmp[, c(columnFields, assayFields)];
        } else {
          tmp <- tmp[, assayFields];
          data <- cbind(data, tmp);
        }
        rm(tmp);
      }
    }

    # If a 'spots' section, write 'count' header to tell how many rows
    # there is the table; this speeds up reading a lot.
    if (headers[["section"]] == "spots") {
      cat(file=con, "count", "\t", nrow(data), "\n", sep="");
    }

    # Write data table separator
    cat(file=con, "%\n");  

    write.table(data, file=con, append=TRUE, sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE);

    cat(file=con, "\n");  
  
    TRUE;
  } # writeBaseFileData()


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
  verbose && enter(verbose, "Writing section");
  on.exit(verbose && exit(verbose), add=TRUE);

  # Each section contains a set of headers and an optional set of
  # data lines, which usually have tab-separated columns.

  writeBaseFileHeaders(con, section);

  writeBaseFileData(con, section);

  invisible(TRUE);
}) # writeBaseFileSection()


############################################################################
# HISTORY: 
# 2005-12-21
# o BUG FIX: When re-reading data cached to file, data cells containing 
#   spaces would generate an error in the internal read.table(). Forgot to
#   specify sep="\t" when reading cached data.  The same was missing 
#   in getData.BaseFileSpots fixed yesterday.
# 2005-07-24
# o Now writeBaseFileSection() adds the header 'count' to 'spots' sections
#   with data tables.
# 2005-07-07
# o BUG FIX: When re-creating data section, the data file names was
#   incorrectly generated on the fly based on the 'assays' header, but
#   this should not be done. Instead, the header 'dataFiles' should be used.
# 2005-06-17
# o Asserts existance for header 'section' (at least).
# o Created.
############################################################################
