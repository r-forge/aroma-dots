readDcpHeader <- function(con, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (is.character(con)) {
    pathname <- con;
    pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);
    con <- file(con, open="rb");
    on.exit({
      if (!is.null(con))
        close(con);
      con <- NULL;
    });
  }

  if (!inherits(con, "connection")) {
    throw("Argument 'con' must be either a connection or a pathname: ", 
                                                            mode(con));
  }

  hdr <- list();

  hdr$Header <- readBin(con=con, what="raw", n=1000);
  hdr$Format <- readBin(con=con, what="integer", size=1, signed=FALSE, n=1);

  hdr$Normalized <- readBin(con=con, what="logical", size=1, n=1);
  hdr$ThetaValid <- readBin(con=con, what="logical", size=1, n=1);

  # Don't know why I have to add this
  hdr$dummy <- readBin(con=con, what="raw", n=1);
# print(seek(con, read="r"));

  hdr$Median <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
  hdr$MaxInten <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
  hdr$CellDim <- readBin(con=con, what="integer", size=4, signed=TRUE, n=1);
# print(seek(con, read="r"));

  hdr$DatFile <- rawToChar(readBin(con=con, what="raw", n=1000));
# print(seek(con, read="r"));
  hdr$BaselineFile <- rawToChar(readBin(con=con, what="raw", n=1000));
# print(seek(con, read="r"));

  hdr$ArrayOutlierPct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
  hdr$SingleOutlierPct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
  hdr$PresencePct <- readBin(con=con, what="double", size=4, signed=TRUE, n=1);
# print(seek(con, read="r"));

  hdr;
} # readDcpHeader()


##############################################################################
# HISTORY:
# 2008-01-30
# o Created.
##############################################################################

