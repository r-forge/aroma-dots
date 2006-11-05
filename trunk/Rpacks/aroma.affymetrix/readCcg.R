#########################################################################/**
# @RdocFunction readCcg
# 
# @title "Reads an Affymetrix Command Console Generic (CCG) Data file"
# 
# @synopsis
# 
# \description{
#   @get "title".  The CCG data file format is also known as the 
#   Calvin file format.
# }
# 
# \arguments{
#   \item{pathname}{The pathname of the CCG file.}
#   \item{...}{Not used.}
#   \item{verbose}{An @integer specifying the verbose level. If 0, the
#     file is parsed quietly.  The higher numbers, the more details.}
# }
# 
# \value{
#   A named @list structure consisting of ...
# }
# 
# @author
# 
#  \seealso{
#    @see "readCdfUnits".
#  }
# 
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
readCcg <- function(pathname, ...) {
  con <- file(pathname, open="rb");
  on.exit(close(con));

  fhdr <- readCcgFileHeader(con);
  dhdr <- readCcgDataHeader(con);

#str(fhdr);
#str(dhdr);

  nextDataGroupStart <- fhdr$dataGroupStart;
  dataGroups <- list();
  for (gg in seq(length=fhdr$nbrOfDataGroups)) {
    dataGroupHeader <- readCcgDataGroupHeader(con, fileOffset=nextDataGroupStart);

    offset <- dataGroupHeader$dataSetStart; 
    dss <- vector("list", dataGroupHeader$nbrOfDataSets);
    names <- character(dataGroupHeader$nbrOfDataSets);
    for (kk in seq(along=dss)) { 
      ds <- readCcgDataSet(con, fileOffset=offset); 
      offset <- ds$nextDataSetStart; 
      dss[[kk]] <- ds;
      names[kk] <- ds$name;
    };
    names(dss) <- names;
  
    dataGroup <- list(
      header = dataGroupHeader,
      dataSets = dss
    );

#str(dataGroup);

    dataGroups <- c(dataGroups, list(dataGroup));
    nextDataGroupStart <- dataGroupHeader$nextGroupStart;
  } # while (nextDataGroupStart != 0)
  names(dataGroups) <- unlist(lapply(dataGroups, FUN=function(dg) dg$header$name), use.names=FALSE);

  ccg <- list(
    fileHeader = fhdr,
    genericDataHeader = dhdr,
    dataGroups = dataGroups
  );
 
  ccg;
} # readCcg()



# File Header
# The file header section is the first section of the file. This 
# section is used to identify the type of file (i.e. Command Console
# data file), its version number (for the file format) and the number
# of data groups stored within the file. Information about the contents
# of the file such as the data type identifier, the parameters used to
# create the file and its parentage is stored within the generic data 
# header section.
# 
# Item 	Description 	Type
# 1 Magic number. A value to identify that this is a Command Console 
#    data file. The value will be fixed to 59. 	[UBYTE]
# 2 The version number of the file. This is the version of the file
#    format. It is currently fixed to 1. [UBYTE]
# 3 The number of data groups. [INT]
# 4 File position of the first data group. [UINT]
readCcgFileHeader <- function(con, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readUByte <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=1, signed=FALSE, endian="big", n=n);
  }

  readInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=TRUE, endian="big", n=n);
  }

  readUInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=FALSE, endian="big", n=n);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  magic <- readUByte(con);
  version <- readUByte(con);
  nbrOfDataGroups <- readInt(con);
  dataGroupStart <- readUInt(con);

  list(
    version = version,
    nbrOfDataGroups = nbrOfDataGroups,
    dataGroupStart = dataGroupStart
  )  
} # readCcgFileHeader()



readCcgDataHeader <- function(con, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readUByte <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=1, signed=FALSE, endian="big", n=n);
  }

  readInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=TRUE, endian="big", n=n);
  }

  readUInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=FALSE, endian="big", n=n);
  }

  readString <- function(con, ...) {
    nchars <- readInt(con);
    if (nchars == 0)
      return("");
    readChar(con, nchars=nchars);
  }

  readWString <- function(con, ...) {
    nchars <- readInt(con);
    if (nchars == 0)
      return("");
    bfr <- readBin(con, what=raw(), n=2*nchars);
    bfr <- bfr[seq(from=2, to=length(bfr), by=2)];
    bfr <- as.integer(bfr);
    bfr <- intToChar(bfr);
    bfr <- paste(bfr, collapse="");
    bfr;
  }

  readWVT <- function(con, ...) {
    list(name=readWString(con), value=readString(con), type=readWString(con));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- list(
    dataTypeId = readString(con),
    fileId = readString(con),
    timestamp = readWString(con),
    locale = readWString(con)
  )

  # Reading parameters
  nbrOfParams <- readInt(con);
  params <- vector("list", nbrOfParams);
  names <- character(nbrOfParams);
  for (kk in seq(length=nbrOfParams)) {
    wvt <- readWVT(con);
    names[kk] <- wvt$name;
    value <- wvt$value;
    attr(value, "mimeType") <- wvt$type;
    params[[kk]] <- value;
  }
  names(params) <- names;
  hdr$parameters <- params;

  # Reading parent headers
  nbrOfParents <- readInt(con);
  parents <- vector("list", nbrOfParents);
  for (kk in seq(length=nbrOfParents)) {
    parents[[kk]] <- readCcgDataHeader(con);
  }
  hdr$parents <- parents;
  
  hdr;
} # readCcgDataHeader()


readCcgDataGroupHeader <- function(con, fileOffset=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=TRUE, endian="big", n=n);
  }

  readUInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=FALSE, endian="big", n=n);
  }

  readWString <- function(con, ...) {
    nchars <- readInt(con);
    if (nchars == 0)
      return("");
    bfr <- readBin(con, what=raw(), n=2*nchars);
    bfr <- bfr[seq(from=2, to=length(bfr), by=2)];
    bfr <- as.integer(bfr);
    bfr <- intToChar(bfr);
    bfr <- paste(bfr, collapse="");
    bfr;
  }

  if (!is.null(fileOffset)) {
    seek(con=con, where=fileOffset, offset="start", rw="read");
  }

  # Data Group
  # This section describes the data group. A data group is a group 
  # of data sets. The file supports one or more data groups in a file.
  # 
  # Item 	Description 	Type
  # 1 	File position of the next data group. When this is the last
  #     data group in the file, the value should be 0. 	UINT
  # 2 	File position of the first data set within the data group. 	UINT
  # 3 	The number of data sets within the data group. 	INT
  # 4 	The data group name. 	WSTRING
    nextGroupStart=readUInt(con)
    dataSetStart=readUInt(con)
    nbrOfDataSets=readInt(con)
    name=readWString(con)

  dataGroupHeader <- list(
    nextGroupStart=nextGroupStart,
    dataSetStart=dataSetStart,
    nbrOfDataSets=nbrOfDataSets,
    name=name
  )

  dataGroupHeader;
} # readCcgDataGroupHeader()




readCcgDataSet <- function(con, fileOffset=NULL, ...) {
  # Value Types
  # The following table defines the numeric values for the value types.
  # The value type is used to representing the type of value stored in 
  # the file.
  #
  # Value 	Type
  # 0 	BYTE
  # 1 	UBYTE
  # 2 	SHORT
  # 3 	USHORT
  # 4 	INT
  # 5 	UINT
  # 6 	FLOAT
  # 7 	DOUBLE
  # 8 	STRING
  # 9 	WSTRING
  whats <- c("integer", "integer", "integer", "integer", "integer", 
            "integer", "double", "double", "character", "character");
  names(whats) <- c("BYTE", "UBYTE", "SHORT", "USHORT", "INT", "UINT", "FLOAT", "DOUBLE", "STRING", "WSTRING");
  signeds <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE);
  sizes <- c(1, 1, 2, 2, 4, 4, 4, 8, 1, 2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readByte <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=1, signed=TRUE, endian="big", n=n);
  }

  readInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=TRUE, endian="big", n=n);
  }

  readUInt <- function(con, n=1, ...) {
    readBin(con, what=integer(), size=4, signed=FALSE, endian="big", n=n);
  }

  readWString <- function(con, ...) {
    nchars <- readInt(con);
    if (nchars == 0)
      return("");
    bfr <- readBin(con, what=raw(), n=2*nchars);
    bfr <- bfr[seq(from=2, to=length(bfr), by=2)];
    bfr <- as.integer(bfr);
    bfr <- intToChar(bfr);
    bfr <- paste(bfr, collapse="");
    bfr;
  }

  readString <- function(con, ...) {
    nchars <- readInt(con);
    if (nchars == 0)
      return("");
    readChar(con, nchars=nchars);
  }

  readWVT <- function(con, ...) {
    list(name=readWString(con), value=readString(con), type=readWString(con));
  }

  readWBI <- function(con, ...) {
    list(name=readWString(con), type=readByte(con), size=readInt(con));
  }


  if (!is.null(fileOffset)) {
    seek(con=con, where=fileOffset, offset="start", rw="read");
  }

#  Data Set
#  This section describes the data for a single data set item
#  (probe set, sequence, allele, etc.). The file supports one
#  or more data sets within a data group.
#  
#  Item 	Description 	Type
#  1 	The file position of the first data element in the data set.
#     This is the first byte after the data set header. 	UINT
#  2 	The file position of the next data set within the data group.
#     When this is the last data set in the data group the value 
#     shall be 1 byte past the end of the data set. This way the size
#     of the data set may be determined. 	UINT
#  3 	The data set name. 	WSTRING
#  4 	The number of name/value/type parameters. 	INT
#  5 	Array of name/value/type parameters. 	(WSTRING / VALUE / TYPE) [ ]
#  6 	Number of columns in the data set.
#     Example: For expression arrays, columns may include signal, p-value,
#     detection call and for genotyping arrays columns may include allele
#     call, and confidence value. For universal arrays, columns may 
#     include probe set intensities and background. 	UINT
#  7 	An array of column names, column value types and column type sizes
#     (one per column).
#     The value type shall be represented by the value from the value type
#     table. The size shall be the size of the type in bytes. For strings,
#     this value shall be the size of the string in bytes plus 4 bytes for
#     the string length written before the string in the file. 	
#     (WSTRING / BYTE / INT) [ ]
#  8 	The number of rows in the data set. 	UINT
#  9 	The data set table, consisting of rows of columns (data values). 
#     The specific type and size of each column is described by the data
#     and size types above. 	ROW [ ]
  dataSet <- list(
    elementsStart=readUInt(con),
    nextDataSetStart=readUInt(con),
    name=readWString(con)
  )

  # Reading parameters
  nbrOfParams <- readInt(con);
  params <- vector("list", nbrOfParams);
  names <- character(nbrOfParams);
  for (kk in seq(length=nbrOfParams)) {
    wvt <- readWVT(con);
    names[kk] <- wvt$name;
    value <- wvt$value;
    attr(value, "mimeType") <- wvt$type;
    params[[kk]] <- value;
  }
  names(params) <- names;
  dataSet$parameters <- params;

  # Reading columns
  nbrOfColumns <- readUInt(con);
  columns <- vector("list", nbrOfColumns);
  names <- character(nbrOfColumns);
  colWhats <- vector("list", nbrOfColumns);
  bytesPerRow <- 0;
  for (cc in seq(length=nbrOfColumns)) {
    wbi <- readWBI(con);
    names[cc] <- wbi$name;
    what <- whats[wbi$type+1];
    signed <- signeds[wbi$type+1];
    size <- wbi$size;
    bytesPerRow <- bytesPerRow + size;
    attr(what, "name") <- names(whats)[wbi$type+1];
    attr(what, "signed") <- signed;
    attr(what, "size") <- size;
    colWhats[[cc]] <- what;
  }
  bytesPerRow <- as.integer(bytesPerRow);

  nbrOfRows <- readUInt(con);
  totalNbrOfBytes <- nbrOfRows * bytesPerRow;

  # Skip to the first element
  seek(con, which=dataSet$elementsStart, offset="start", rw="read");
  # Read all data row by row
  raw <- readBin(con, what=raw(), n=totalNbrOfBytes);
  dim(raw) <- c(bytesPerRow, nbrOfRows);

  table <- vector("list", nbrOfColumns);
  colsOffset <- 0;
  for (cc in seq(length=nbrOfColumns)) {
    what <- colWhats[[cc]];
    signed <- attr(what, "signed");
    size <- attr(what, "size");

    if (what == "character") {
      rawCol <- matrix(raw[1:4,], nrow=nbrOfRows, ncol=4);
      nchars <- readInt(con=rawCol, n=nbrOfRows);
      nchars <- nchars[1];
      ccs <- colsOffset + 4 + seq(from=2, to=2*nchars, by=2);
      value <- intToChar(as.integer(raw[ccs,]));
      dim(value) <- c(nchars, nbrOfRows);
      value <- apply(value, MARGIN=2, paste, collapse="");
    } else {
      ccs <- colsOffset + 1:size;
      rawCol <- matrix(raw[ccs,], nrow=bytesPerRow, ncol=nbrOfRows, byrow=FALSE);
      value <- readBin(con=rawCol, what=what, size=size, signed=signed, endian="big", n=nbrOfRows);
    }
    table[[cc]] <- value;
    colsOffset <- colsOffset + size;
  }
  table <- as.data.frame(table);
  colnames(table) <- names;
  dataSet$table <- table;

  dataSet;
} # readCcgDataSet()


############################################################################
# HISTORY:
# 2006-11-06
# o Tested on Test3-1-121502.calvin.CEL and Test3-1-121502.calvin.CDF.
# o Created.  
############################################################################  
