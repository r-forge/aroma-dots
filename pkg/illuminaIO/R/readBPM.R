readBPM <- function(filename, path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readByte <- function(con, n=1, ...) {
    readBin(con, what="integer", n=n, size=1, endian="little", signed=FALSE);
  }

  readShort <- function(con, n=1, ...) {
    readBin(con, what="integer", n=n, size=2, endian="little", signed=FALSE);
  }

  readInt <- function(con, n=1, ...) {
    readBin(con, what="integer", n=n, size=4, endian="little", signed=TRUE);
  }

  readLong <- function(con, n=1, ...) {
    readBin(con, what="integer", n=n, size=8, endian="little", signed=TRUE);
  }

  readString <- function(con, ...) {
    n <- readByte(con, n=1);
    readChar(con, nchars=n);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  fileSize <- file.info(pathname)$size;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Open file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  con <- file(pathname, "rb");
  on.exit({
    close(con);
  });

  ## Reads and parses Illumina BPM files

  # The first few bytes of the egtFile are some type of
  # header, but there's no related byte offset information.

  # Assert file format
  magic <- readChar(con, nchars=3);
  if (magic != "BPM") {
    throw("Cannot read BPM file. File format error. Unknown magic: ", magic);
  }

  null.1 <- readByte(con, n=1);
  ## should be 1

  # Read file format version
  version <- readLong(con, n=1);
  if (version != 4) {
    throw("Cannot read BPM file. Unsupported IDAT file format version: ", version);
  }

  chipType <- readString(con);

  null.2 <- readByte(con, n=2);

  csvLines <- readLines(con, n=22);

  entriesByteOffset <- seek(con);
  nEntries <- readInt(con, n=1);

  if (FALSE) {
    snpIndexByteOffset <- seek(con);
    snpIndex <- readInt(con, n=nEntries);
    ## for the 1M array, these are simply in order from 1 to 1072820.

    snpNamesByteOffset <- seek(con);
    snpNames <- rep("A", times=nEntries);
    for (ii in seq(length=nEntries)) {
      snpNames[ii] <- readString(con);
    }
  }

  seek(con, where=15278138, origin="start");

  normIDByteOffset <- seek(con);
  normID <- readByte(con, n=nEntries) + 1;

  newBlockByteOffset <- seek(con);
  newBlock <- readByte(con, n=10000);

  byteOffsets <- list(
    entriesByteOffset=entriesByteOffset,
    #snpIndexByteOffset=snpIndexByteOffset,
    #snpNamesByteOffset=snpNamesByteOffset,
    normIDByteOffset=normIDByteOffset,
    newBlockByteOffset=newBlockByteOffset
  );

  res <- list(
    prefixCheck=magic,
    null.1=null.1,
    versionNumber=version,
    chipType=chipType,
    null.2=null.2,
    csvLines=csvLines,
    nEntries=nEntries,
    #snpIndex=snpIndex,
    #snpNames=snpNames,
    normID=normID,
    newBlock=newBlock,
    byteOffsets=byteOffsets
  );

  res
} # readBPM()


############################################################################
# HISTORY:
# 2011-03-28
# o CLEAN UP: Tidied up the code.
# o ROBUSTNESS: Now file connection is closed via on.exit().
# o Created from crlmm/R/crlmm-illumina.R, which says that this
#   function was original written and provided by Keith Baggerly
#   on 2008-08-27.
############################################################################ 

