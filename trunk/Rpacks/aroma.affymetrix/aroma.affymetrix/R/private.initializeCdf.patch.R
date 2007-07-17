.initializeCdf.patch <- function(con, nRows = 1, nCols = 1,
                          nUnits = 1, nQcUnits = 0,
                          refSeq = "",
                          unitnames = rep("", nUnits),
                          qcUnitLengths = rep(0, nQcUnits),
                          unitLengths = rep(0, nUnits),
                          ...) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(length(qcUnitLengths) != nQcUnits) {
      stop("Number of elements in argument 'qcUnitLengths' does not match 'nQcUnits'");
    }

    if(length(unitLengths) != nUnits) {
      stop("Number of elements in argument 'qcUnitLengths' does not match 'nUnits'");
    }

    if(length(refSeq) != 1)
        stop("Argument 'refSeq' should be a single character.");

    lrefSeq <- nchar(refSeq);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # CDF header
    #
    # 1 Magic number. Always set to 67.                           [integer]
    # 2 Version number.                                           [integer]
    # 3 The number of columns of cells on the array.       [unsigned short]
    # 4 The number of rows of cells on the array.          [unsigned short]
    # 5 The number of units in the array not including QC units. The term 
    #   unit is an internal term which means probe set.           [integer]
    # 6 The number of QC units.                                   [integer]
    # 7 The length of the resequencing reference sequence.        [integer]
    # 8 The resequencing reference sequence.                    [char[len]]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    offset <- 0;
    
    ## Magic number and version number
    writeBin(object = as.integer(c(67, 1)),
             con = con, size = 4, endian = "little")
    ## NCols, NRows
    writeBin(object = as.integer(c(nCols, nRows)),
             con = con, size = 2, endian = "little")
    ## NumberUnits, NumberQCUnits
    writeBin(object = as.integer(c(nUnits, nQcUnits)),
             con = con, size = 4, endian = "little")
    ## Length of refSeqsequence
    writeBin(object = as.integer(lrefSeq),
             con = con, size = 4, endian = "little")
    offset <- 24;

    fOffset <- seek(con=con, origin="start", rw="write");
    if (offset != fOffset) {
      throw("File format write error (step 1): File offset is not the excepted one: ", fOffset, " != ", offset);
    }   
 
    ## RefSeqsequece
    if(lrefSeq > 0)
      writeChar(as.character(refSeq), con=con, eos=NULL);

    # Current offset
    offset <- offset + lrefSeq;

    fOffset <- seek(con=con, origin="start", rw="write");
    if (offset != fOffset) {
      throw("File format write error (step 2): File offset is not the excepted one: ", fOffset, " != ", offset);
    }   


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Unit names
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Write to raw vector (2*10^6 units => 122Mb; should be ok for now)
    # Since we can't create strings with '\0':s, we use '\xFF',
    # write to raw and then replace '\xFF' with '\0'. Thus, unit names with
    # '\xFF' are invalid, but this should not be a real problem.
    pads <- sapply(0:64, FUN=function(x) paste(rep("\xFF", x), collapse=""));

    # Write the unit names in chunks to save memory
    nbrOfUnits <- length(unitnames);
    chunkSize <- 100000;
    nbrOfChunks <- ceiling(nbrOfUnits / chunkSize);

    # Allocate raw vector
    raw <- raw(64*chunkSize);

    nunits2 <- 0;
    for (kk in 1:nbrOfChunks) {
      # Units for this chunk
      from <- (kk-1)*chunkSize+1;
      to <- min(from+chunkSize-1, nbrOfUnits);
      idxs <- from:to;
      nunits2 <- nunits2 + length(idxs);
      unitnamesFF <- unitnames[idxs];

      # Pad the unit names
      unitnamesFF <- paste(unitnamesFF, pads[64-nchar(unitnamesFF)], sep="");

      # Truncate last chunk?
      if (chunkSize > length(unitnamesFF)) {
        raw <- raw[1:(64*length(unitnamesFF))];
      }

      # Write unit names to raw vector
      raw <- writeBin(con=raw, unitnamesFF, size=1);

      rm(unitnamesFF);

      # Garbage collect
#      gc <- gc();
#      print(gc);

      # Replace all '\xFF' with '\0'.
      idxs <- which(raw == as.raw(255));
      raw[idxs] <- as.raw(0);
      rm(idxs);

#    fOffset <- seek(con=con, origin="start", rw="write");
#cat(sprintf("%d. length(raw)=%d, fOffset=%d\n", kk, length(raw), fOffset));

      writeBin(con=con, raw);
   } # for (kk in ...)

#    fOffset <- seek(con=con, origin="start", rw="write");
#cat(sprintf("ZZZ. fOffset=%d\n", fOffset));

   rm(raw);
   # Garbage collect
   gc <- gc();

#    writeChar(con=con, as.character(unitnames), nchars=rep(64, nUnits), eos=NULL)

    bytesOfUnitNames <- 64 * nUnits;
    offset <- offset + bytesOfUnitNames;

    fOffset <- seek(con=con, origin="start", rw="write");
    if (offset != fOffset) {
      throw("File format write error (step 3): File offset is not the excepted one: ", fOffset, " != ", offset);
    }   

    bytesOfQcUnits <- 4 * nQcUnits;
    offset <- offset + bytesOfQcUnits;

    bytesOfUnits <- 4 * nUnits;
    offset <- offset + bytesOfUnits;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # QC units file positions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nQcUnits > 0) {
      csum <- cumsum(qcUnitLengths);
      nextOffset <- csum[nQcUnits];
      starts <- c(0, csum[-nQcUnits]);
      starts <- as.integer(offset + starts);
      writeBin(starts, con = con, size = 4, endian = "little")
    } else {
      nextOffset <- 0;
#      starts <- 0;
#      starts <- as.integer(offset + starts);
#      writeBin(starts, con = con, size = 4, endian = "little")
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Units file positions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    offset <- offset + nextOffset;
    if (nUnits > 0) {
      csum <- cumsum(unitLengths);
      nextOffset <- csum[nUnits];
      starts <- c(0, csum[-nUnits]);
      starts <- as.integer(offset + starts);
      writeBin(starts, con = con, size = 4, endian = "little");
    } else {
      nextOffset <- 0;
    }
} # .initializeCdf()
