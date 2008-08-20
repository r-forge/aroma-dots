###########################################################################/**
# @RdocFunction readCdfBinHeader
#
# @title "Reads the file header of a dChip CDF.bin file"
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
#   To read the CDF.bin file data, see @see "readCdfBin".
# }
#
# @keyword "file"
# @keyword "IO"
#*/########################################################################### 
readCdfBinHeader <- function(con, ...) {
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

  hdr <- list();

  LINE_LEN <- 1000;
  UNIT_NAME_LEN <- 50;
  hdr$FileName <- paste(rawToChar(readBin(con=con, what="raw", n=LINE_LEN-1)));
  hdr$Format <- readBin(con=con, what="integer", size=1, signed=FALSE, n=1);
  hdr$ChipType <- paste(rawToChar(readBin(con=con, what="raw", n=UNIT_NAME_LEN)));

  # Don't know why I have to add this
  hdr$dummy <- readBin(con=con, what="raw", n=2);
	
  hdr$CellDim <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
  hdr$NumUnit <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up the header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr$dummy <- NULL; 

  hdr;
} # readCdfBinHeader()



##############################################################################
# HISTORY:
# 2008-08-20
# o Added Rdoc comments.
# 2008-02-03
# o Created.
##############################################################################
