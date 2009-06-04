###########################################################################/**
# @RdocFunction readDcpHeader
#
# @title "Reads the file header of a dChip DCP file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character filename.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list structure containing the file header.
# }
#
# @author
#
# \seealso{
#   To read also the DCP file data, see @see "readDcp".
# }
#
# @keyword "file"
# @keyword "IO"
#*/########################################################################### 
readDcpHeader <- function(con, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rawToString <- function(raw, ...) {
    # This approach drops all '\0', in order to avoid warnings
    # in rawToChar().  Note, it does not truncate the string after
    # the first '\0'.  However, such strings should never occur in
    # the first place.
    raw <- raw[raw != as.raw(0)];
    rawToChar(raw);
  } # rawToString()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (is.character(con)) {
    pathname <- con;
    if (!file.exists(pathname)) {
      stop("File not found: ", pathname);
    }
    con <- file(con, open="rb");
    on.exit({
      if (!is.null(con))
        close(con);
      con <- NULL;
    });
  }

  if (!inherits(con, "connection")) {
    stop("Argument 'con' must be either a connection or a pathname: ", 
                                                            mode(con));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read the file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- list();

  hdr$Header <- readBin(con=con, what="raw", n=1000);
  hdr$Format <- readBin(con=con, what="integer", size=1, signed=FALSE, n=1);

  # Assert that the file format is correct/the expected on
  if (!is.element(hdr$Format, 3:4)) {
    stop("File format error: The DCP header format ('Format') is not v3 or v4: ", hdr$Format);
  }

  hdr$Normalized <- readBin(con=con, what="logical", size=1, n=1);
  hdr$ThetaValid <- readBin(con=con, what="logical", size=1, n=1);

  # Don't know why we have to add this one
  hdr$dummy <- readBin(con=con, what="raw", n=1);
# print(seek(con, read="r"));

  hdr$Median <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
  hdr$MaxInten <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
  hdr$CellDim <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
# print(seek(con, read="r"));

  hdr$DatFile <- rawToString(readBin(con=con, what="raw", n=1000));
# print(seek(con, read="r"));
  hdr$BaselineFile <- rawToString(readBin(con=con, what="raw", n=1000));
# print(seek(con, read="r"));

  hdr$ArrayOutlierPct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
  hdr$SingleOutlierPct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
  hdr$PresencePct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
# print(seek(con, read="r"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up the header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr$dummy <- NULL; 
  hdr$Header <- paste(rawToString(hdr$Header));
  hdr$DatFile <- paste(hdr$DatFile);
  hdr$BaselineFile <- paste(hdr$BaselineFile);

  hdr;
} # readDcpHeader()


##############################################################################
# HISTORY:
# 2009-02-13
# o Using rawToString() instead of rawToChar() to avoid warnings on
#   'truncating string with embedded nul:...'.
# 2008-xx-xx
# o Added Rdoc comments.
# o Added argument 'clean'.
# 2008-01-30
# o Created.
##############################################################################

