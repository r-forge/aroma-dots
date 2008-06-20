setConstructorS3("AromaProbeSequenceTextFile", function(...) {
  extend(AromaMicroarrayDataFile(...), "AromaProbeSequenceTextFile");
})


setMethodS3("allocate", "AromaProbeSequenceTextFile", function(this, filename, path=NULL, platform, chipType, nbrOfRows, nbrOfColumns, probeLengths, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'nbrOfRows':
  nbrOfRows <- Arguments$getInteger(nbrOfRows, range=c(1, 1e6));

  # Argument 'nbrOfColumns':
  nbrOfColumns <- Arguments$getInteger(nbrOfColumns, range=c(1, 1e6));

  # Argument 'probeLengths':
  probeLengths <- Arguments$getInteger(probeLengths, range=c(1, 10e3));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # File header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Generating file header");

  nbrOfCells <- nbrOfRows*nbrOfColumns;
  probeLengths <- 25;
  probeSep <- "0x10";

  pkg <- Package("aroma.affymetrix");
  createdBy <- sprintf("%s v%s", getName(pkg), getVersion(pkg));
  createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
  hdr <- list(
    platform = platform,
    chipType = chipType,
    nbrOfCells = nbrOfCells, 
    probeLengths = probeLengths,
    probeSep = probeSep,
    nbrOfRows = nbrOfRows,
    nbrOfColumns = nbrOfColumns,
    createdBy = createdBy,
    createdOn = createdOn
  );

  args <- list(...);
  for (key in names(args)) {
    hdr[[key]] <- args[[key]];
  }

  hdr <- paste("# ", names(hdr), ": ", sapply(hdr, FUN=paste, collapse="\t"), sep="");
  hdr <- paste(hdr, collapse="\n");
  hdr <- paste(hdr, "\n", sep="");
  verbose && cat(verbose, hdr, newline=FALSE);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Write to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Write as a binary file to avoid newline convertions
  con <- file(pathname, open="wb");
  on.exit({
    close(con);
  });

  verbose && enter(verbose, "Allocating probe sequence file");
  # Empty sequence string with newline ending
  emptySeq <- paste(c(rep(" ", probeLengths), "\n"), collapse="");
  emptySeq <- charToRaw(emptySeq);

  # Write header
  hdr <- charToRaw(hdr);
  writeBin(con=con, hdr);

  # Write empty strings in chunks
  CHUNK.SIZE <- 500e3;
  emptyChunk <- rep(emptySeq, times=CHUNK.SIZE);
  while (nbrOfCells > 0) {
    if (nbrOfCells < CHUNK.SIZE) {
      emptyChunk <- rep(emptySeq, times=nbrOfCells);
    }
    writeBin(con=con, emptyChunk);
    nbrOfCells <- nbrOfCells - CHUNK.SIZE;
  }

  verbose && exit(verbose);
  
  res <- AromaProbeSequenceTextFile(pathname);
  res;
})


setMethodS3("findByChipType", "AromaProbeSequenceTextFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);


  fullname <- paste(c(chipType, tags), collapse = ",");

  pattern <- sprintf("^%s.*,probeSeqs.(txt|TXT)$", fullname);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;
  pathname <- do.call("findAnnotationDataByChipType", args=args);
  if (is.null(pathname)) {
    pattern <- sprintf("^%s.*,probeSeqs.(txt|TXT)[.]lnk$", fullname);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE)


setMethodS3("byChipType", "AromaProbeSequenceTextFile", function(static, chipType, tags=NULL, ...) {
  pathname <- findByChipType(static, chipType=chipType, tags=tags, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix tabular file: ", chipType);

  newInstance(static, pathname, ...);
}, static=TRUE)



setMethodS3("cleanupNewlines", "AromaProbeSequenceTextFile", function(this, ...) {
  pathname <- getPathname(this);
  destPathname <- sprintf("%s.tmp", pathname);
  if (isFile(destPathname)) {
#    throw("Cannot cleanup newlines. Temporary file already exists: ", destPathname);
  }

  con <- destCon <- NULL;
  on.exit({
    if (!is.null(con)) {
      close(con);
      con <- NULL;
    }
    if (!is.null(destCon)) {
      close(destCon);
      destCon <- NULL;
    }
  })

  # Open file connections
  con <- file(pathname, open="rb");
  destCon <- file(destPathname, open="wb");

  carriageReturn <- charToRaw("\r");
  newline <- charToRaw("\n");
  CHUNK.SIZE <- 1024^2;
  while(TRUE) {
    bfr <- readBin(con=con, what="raw", n=CHUNK.SIZE);
    if (length(bfr) == 0)
      break;
    # Identify the positions of all '\r' followed by a '\n'
    pos <- which(bfr == carriageReturn);
    if (length(pos) > 0) {
      keep <- (bfr[pos+1] == newline);
      pos <- pos[keep];
      # Exclude those
      bfr <- bfr[-pos];
    }
    # Write
    writeBin(con=destCon, bfr);
  }

  close(con); con <- NULL;
  close(destCon); destCon <- NULL;

  if (!file.rename(destPathname, pathname)) {
    throw("Failed to rename temporary cleaned up file: ", 
                                    destPathname, " -> ", pathname);
  }

  invisible(this);
}, protected=TRUE);


setMethodS3("readTableHeader", "AromaProbeSequenceTextFile", function(this, con=NULL, ...) {

  # Open a file connection?
  if (is.null(con)) {
    pathname <- getPathname(this);

    # Open file connection
    con <- file(pathname, open="r");
    on.exit({
      if (!is.null(con)) {
        close(con);
        con <- NULL;
      }
    })
  }

  hdr <- readTableHeader(con=con, ...);
  hdr$topRows <- lapply(hdr$topRows, trim);

  hdr;
}) # readTableHeader()


setMethodS3("getDataFileOffset", "AromaProbeSequenceTextFile", function(this, ...) {
  pathname <- getPathname(this);
  con <- file(pathname, open="rb");
  on.exit({
    if (!is.null(con)) {
      close(con);
      con <- NULL;
    }
  })

  commentChar <- "#";
  newlineChars <- c("\n", "\r");
  rawCommentChar <- charToRaw(commentChar);
  rawNewlineChars <- charToRaw(paste(newlineChars, collapse=""));
  # 1. Scan for the last comment row
  state <- "start";
  pos <- 0;
  while(TRUE) {
    ch <- readBin(con=con, what="raw", n=1);
    if (state == "start") {
      if (ch %in% rawNewlineChars) {
      } else if (ch == rawCommentChar) {
        state <- "inComment";
      } else {
        break;
      }
    } else if (state == "inComment") {
      if (ch %in% rawNewlineChars)
        state <- "start";
    }
    pos <- pos + 1;
  } # while()

  pos <- as.integer(pos);
  pos;
}, protected=TRUE)


setMethodS3("getHeaderParameters", "AromaProbeSequenceTextFile", function(this, ...) {
  hdr <- readTableHeader(this);
  pattern <- "#[ ]*([^:= ]*)[ ]*[:=][ ]*(.*)";
  keys <- gsub(pattern, "\\1", hdr$comments);
  keys <- trim(keys);
  values <- gsub(pattern, "\\2", hdr$comments);

  params <- strsplit(values, split="\t");
  names(params) <- keys;

  params$nbrOfCells <- as.integer(params$nbrOfCells);
  params$probeLengths <- as.integer(params$probeLengths);
  params$nbrOfColumns <- as.integer(params$nbrOfColumns);
  params$nbrOfRows <- as.integer(params$nbrOfRows);

  params;
})

setMethodS3("getPlatform", "AromaProbeSequenceTextFile", function(this, ...) {
  params <- getHeaderParameters(this, ...);
  params$platform;
})

setMethodS3("getChipType", "AromaProbeSequenceTextFile", function(this, ...) {
  params <- getHeaderParameters(this, ...);
  params$chipType;
})

setMethodS3("nbrOfCells", "AromaProbeSequenceTextFile", function(this, ...) {
  params <- getHeaderParameters(this, ...);
  params$nbrOfCells;
})

setMethodS3("nbrOfColumns", "AromaProbeSequenceTextFile", function(this, ...) {
  params <- getHeaderParameters(this, ...);
  params$nbrOfColumns;
})

setMethodS3("nbrOfRows", "AromaProbeSequenceTextFile", function(this, ...) {
  params <- getHeaderParameters(this, ...);
  params$nbrOfRows;
})


setMethodS3("getCellsFromXY", "AromaProbeSequenceTextFile", function(this, x, y=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'x':
  if (is.matrix(x) || is.data.frame(x)) {
    y <- as.integer(x[,"y"]);
    x <- as.integer(x[,"x"]);
  }

  nbrOfColumns <- nbrOfColumns(this);
  cells <- nbrOfColumns*y + x + 1;
  cells <- as.integer(cells);
  cells;
}, protected=TRUE)



setMethodS3("isMissing", "AromaProbeSequenceTextFile", function(this, cells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(cells)) {
  } else {
    if (is.matrix(cells) || is.data.frame(cells)) {
      cells <- getCellsFromXY(this, cells);
    }
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
    nbrOfCells <- length(cells);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Compatibility check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathname <- getPathname(this);
  params <- getHeaderParameters(this);
  probeLengths <- params$probeLengths;
  if (is.na(probeLengths)) {
    throw("Cannot read probe sequences: Variable length probe sequences is not supported: ", pathname);
  }


  con <- file(pathname, open="rb");
  on.exit({
    close(con);
  })

  # Move to start of data section
  dataOffset <- getDataFileOffset(this);
  seek(con=con, where=dataOffset, origin="start", rw="read");

  if (is.null(cells)) {
    rCells <- c(1, nbrOfCells);
  } else {
    rCells <- range(cells);
  }

  allFirst <- raw(nbrOfCells);
  nbrOfCellsRead <- 100e3;
  seqLength <- probeLengths+1;
  CHUNK.SIZE <- seqLength*nbrOfCellsRead;
  firstBase <- seq(from=1, to=CHUNK.SIZE, by=seqLength);
  cellOffset <- nextCell <- 0;
  while (TRUE) {
    # No more cells to read?
    if (nextCell > rCells[2])
      break;

    # Read big chunk of sequence data
    bfr <- readBin(con=con, what="raw", n=CHUNK.SIZE);
    n <- length(bfr);
    if (n == 0)
      break;

    # Keep only the first base in each sequence
    if (n < CHUNK.SIZE) {
      firstBase <- seq(from=1, to=n, by=seqLength);
      nbrOfCellsRead <- length(firstBase);
      idxs <- 1:nbrOfCellsRead;
    }
    bfr <- bfr[firstBase];

    # Keep only the sequences asked for
    idxs <- 1:nbrOfCellsRead;
    if (is.null(cells)) {
    } else {
      keep <- match(nextCell + idxs, cells);
      keep <- which(is.finite(keep));
      idxs <- idxs[keep];
      bfr <- bfr[keep];
    }

    allFirst[cellOffset+idxs] <- bfr;

    cellOffset <- cellOffset + length(idxs);
    nextCell <- nextCell + nbrOfCellsRead;
  } # while()
  
  (allFirst == charToRaw(" "));
})


# \item{cells}{Non-duplicated ordered cell indices.}
setMethodS3("readRawSequences", "AromaProbeSequenceTextFile", function(this, cells=NULL, ..., verbose=FALSE) {
  params <- getHeaderParameters(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(cells)) {
  } else {
    if (is.matrix(cells) || is.data.frame(cells)) {
      cells <- getCellsFromXY(this, cells);
    }
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
    nbrOfCells <- length(cells);
  }
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Compatibility check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathname <- getPathname(this);
  probeLengths <- params$probeLengths;
  if (is.na(probeLengths)) {
    throw("Cannot read probe sequences: Variable length probe sequences is not supported: ", pathname);
  }

  verbose && enter(verbose, "Reading sequences");
  verbose && cat(verbose, "Pathname: ", pathname);

  # Identify beginning of data section
  dataOffset <- getDataFileOffset(this);
  verbose && cat(verbose, "Data section offset: ", dataOffset);

  verbose && cat(verbose, "Cell indices:");
  verbose && str(verbose, cells);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify byte intervals to be read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify cell-index intervals
  if (is.null(cells)) {
    fromTo <- matrix(c(1,nbrOfCells), ncol=2);
  } else {
    fromTo <- seqToIntervals(cells);
  }

  # Translate to intervals in bytes
  seqLength <- probeLengths+1;
  fromTo <- seqLength*(fromTo-1)+1;
  fromTo[,2] <- fromTo[,2] + probeLengths;
  fromTo <- dataOffset + fromTo;
  verbose && cat(verbose, "Number of intervals: ", nrow(fromTo));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read raw data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  bfr <- readBinFragments(getPathname(this), what="raw", idxs=fromTo);

  verbose && exit(verbose);

  bfr;
})


setMethodS3("updateRawSequences", "AromaProbeSequenceTextFile", function(this, cells=NULL, seqs, ..., verbose=FALSE) {
  params <- getHeaderParameters(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  probeLengths <- params$probeLengths;
  seqLength <- probeLengths+1;
  if (is.null(cells)) {
  } else {
    if (is.matrix(cells) || is.data.frame(cells)) {
      cells <- getCellsFromXY(this, cells);
    }
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
    nbrOfCells <- length(cells);
  }

  # Argument 'seqs':
  seqs <- Arguments$getVector(seqs, length=rep(seqLength*nbrOfCells, 2));
  seqs <- as.vector(seqs);
  knownValues <- charToRaw("ACGT \n");
  if (!all(seqs %in% knownValues)) {
    unknown <- seqs[!seqs %in% knownValues];
    unknown <- unknown[seq(length=min(length(unknown),5))];
    unknown <- sprintf("%#02x (%s)", as.integer(unknown), rawToChar(unknown, multiple=TRUE));
    unknown <- paste(unknown, collapse=", ");
    throw("Argument 'seqs' contains values not in [ACGT \n]: ", unknown);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Compatibility check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathname <- getPathname(this);
  if (is.na(probeLengths)) {
    throw("Cannot read probe sequences: Variable length probe sequences is not supported: ", pathname);
  }

  verbose && enter(verbose, "Updating sequences");
  verbose && cat(verbose, "Pathname: ", getPathname(this));
    
  # Identify beginning of data section
  dataOffset <- getDataFileOffset(this);
  verbose && cat(verbose, "Data section offset: ", dataOffset);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Ordering objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(cells)) {
    verbose && enter(verbose, "Ordering data");
    o <- order(cells);
    cells <- cells[o];
    seqs <- matrix(seqs, nrow=seqLength);    
    seqs <- seqs[,o];
    seqs <- as.vector(seqs);
    rm(o);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  verbose && cat(verbose, "Sequences:");
  verbose && str(verbose, seqs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify byte intervals to be read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify cell-index intervals
  if (is.null(cells)) {
    fromTo <- matrix(c(1,nbrOfCells), ncol=2);
  } else {
    fromTo <- seqToIntervals(cells);
  }
  verbose && str(verbose, fromTo);

  # Translate to intervals in bytes
  fromTo <- seqLength*(fromTo-1)+1;
  fromTo[,2] <- fromTo[,2] + probeLengths;
  fromTo <- dataOffset + fromTo;
  verbose && cat(verbose, "Number of intervals: ", nrow(fromTo));
  verbose && str(verbose, fromTo);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read raw data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  writeBinFragments(con=getPathname(this), seqs, idxs=fromTo);

  verbose && exit(verbose);
})


setMethodS3("readRawSequenceMatrix", "AromaProbeSequenceTextFile", function(this, cells=NULL, ...) {
  params <- getHeaderParameters(this);

  # Reorder cells?
  if (!is.null(cells)) {
    o <- order(cells);
  }

  seqs <- readRawSequences(this, cells=cells, ...);

  by <- params$probeLengths+1;
  excl <- seq(from=by, to=length(seqs), by=by);
  seqs <- seqs[-excl];
  rm(excl);
  seqs <- matrix(seqs, ncol=params$probeLengths, byrow=TRUE);

  # Re-reorder sequences?
  if (!is.null(cells)) {
    o <- order(o);
    seqs <- seqs[o,,drop=FALSE];
    rm(o);
  }

  seqs;
}, protected=TRUE)



setMethodS3("seqRawMatrixToStrings", "matrix", function(seqs, ...) {
  seqsStrs <- NULL;
  for (cc in seq(length=ncol(seqs))) {
    seqsCc <- rawToChar(seqs[,cc], multiple=TRUE);
    seqsStrs <- paste(seqsStrs, seqsCc, sep="");
  }
  probeLengths <- nchar(seqsStrs[1]);
  unknownSeq <- paste(rep(" ", probeLengths), collapse="");
  seqsStrs[seqsStrs == unknownSeq] <- NA;
  seqsStrs;
}, protected=TRUE)


setMethodS3("seqMatrixToStrings", "matrix", function(seqs, ...) {
  seqsStrs <- NULL;
  for (cc in seq(length=ncol(seqs))) {
    seqsStrs <- paste(seqsStrs, seqs[,cc], sep="");
  }
  probeLengths <- nchar(seqsStrs[1]);
  unknownSeq <- paste(rep(" ", probeLengths), collapse="");
  seqsStrs[seqsStrs == unknownSeq] <- NA;
  seqsStrs;
}, protected=TRUE)


setMethodS3("readSequenceStrings", "AromaProbeSequenceTextFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Read raw probe sequence matrix");
  seqsM <- readRawSequenceMatrix(this, ..., verbose=verbose);
  verbose && exit(verbose);

  verbose && enter(verbose, "Converting to strings");
  seqs <- seqRawMatrixToStrings(seqsM, ...);
  verbose && exit(verbose);

  seqs;
})

setMethodS3("updateSequenceStrings", "AromaProbeSequenceTextFile", function(this, ..., seqs) {
  probeLenghts <- nchar(seqs);
  uProbeLenghts <- unique(probeLenghts);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compatibility checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(uProbeLenghts) != 1) {
    throw("Cannot write probe sequences. Sequences of varying lengths detected.");
  }

  params <- getHeaderParameters(this);
  probeLengths <- params$probeLengths;
  if (uProbeLenghts != probeLengths) {
    throw("Cannot write probe sequences. The sequences to be written are of different lengths than the ones on file: ", uProbeLenghts, " != ", probeLengths);
  }

  # Translate NAs to unknown sequences
  unknownSeq <- paste(rep(" ", probeLengths), collapse="");
  seqs[is.na(seqs)] <- unknownSeq;

  seqs <- paste(seqs, collapse="\n");
  seqs <- charToRaw(seqs);
  seqs <- c(seqs, charToRaw("\n"));

  updateRawSequences(this, ..., seqs=seqs);
})



setMethodS3("countBases", "AromaProbeSequenceTextFile", function(this, cells=NULL, ..., verbose=FALSE) {
  require("matchprobes") || throw("Package not loaded: matchprobes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(cells)) {
  } else {
    if (is.matrix(cells) || is.data.frame(cells)) {
      cells <- getCellsFromXY(this, cells);
    }
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
    nbrOfCells <- length(cells);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Optimize reading order?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(cells)) {
    o <- order(cells);
    cells <- cells[o];
    reorder <- TRUE;
    gc <- gc();
    verbose && print(verbose, gc);
  } else {
    cells <- 1:nbrOfCells;
    reorder <- FALSE;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Count in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  counts <- matrix(NA, nrow=nbrOfCells, ncol=4);
  colnames(counts) <- colnames(matchprobes::countbases("A"));

  CHUNK.SIZE <- 1e6;
  offset <- 0;
  chunk <- 1;
  nbrOfChunks <- ceiling(length(cells) / CHUNK.SIZE);
  while (length(cells) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    # Identify cells to read in this chunk
    n <- min(CHUNK.SIZE, length(cells));
    idxs <- 1:n;
    cellsChunk <- cells[idxs];

    verbose && cat(verbose, "Cells to be read:");
    verbose && str(verbose, cellsChunk);

    # Read sequences
    seqs <- readSequenceStrings(this, cells=cellsChunk, ..., 
                                                verbose=less(verbose, 25));
    rm(cellsChunk);
    verbose && cat(verbose, "Sequences read:");
    verbose && str(verbose, seqs);

    # Read fewer cells?
    n2 <- length(seqs);
    if (n2 < n) {
      n <- n2;
      idxs <- 1:n;
    }

    # Count bases
    countsChunk <- matchprobes::countbases(seqs);
    rm(seqs);
    verbose && cat(verbose, "Counts:");
    verbose && str(verbose, countsChunk);

    # Store
    counts[offset+idxs,] <- countsChunk;

    # Next set of cells
    offset <- offset + n;
    cells <- cells[-idxs];
    rm(idxs);

    gc <- gc();
    verbose && print(verbose, gc);

    chunk <- chunk + 1;
    verbose && exit(verbose);
  } # while()
  rm(cells);

  # Reorder?
  if (reorder) {
    o <- order(o);
    counts <- counts[o,,drop=FALSE];
    rm(o);
    gc <- gc();
    verbose && print(verbose, gc);
  }

  counts;
}, protected=TRUE)


setMethodS3("readSequenceMatrix", "AromaProbeSequenceTextFile", function(this, ...) {
  seqs <- readRawSequenceMatrix(this, ...);
  dim <- dim(seqs);
  seqs <- rawToChar(seqs, multiple=TRUE);
  seqs[seqs == " "] <- as.character(NA);
  dim(seqs) <- dim;
  seqs;
})




setMethodS3("importFrom", "AromaProbeSequenceTextFile", function(this, srcFile, ...) {
  methodName <- sprintf("importFrom%s", class(srcFile)[1]);
  fcn <- get(methodName, mode="function");
  fcn(this, srcFile=srcFile, ...);
})


setMethodS3("importFromAffymetrixProbeTabFile", "AromaProbeSequenceTextFile", function(this, srcFile, rows=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'srcFile':
  if (inherits(srcFile, "AffymetrixProbeTabFile")) {
    # Validate chip type
    chipType <- getChipType(this, fullname=FALSE);
    chipTypeSrc <- getChipType(srcFile, fullname=FALSE);
    if (!identical(chipTypeSrc, chipType)) {
      throw("Argument 'srcFile' has a different chip type: ", 
                                                    chipTypeSrc, chipType);
    }
  } else {
    throw("Argument 'srcFile' is not an AffymetrixProbeTabFile: ", 
                                                        class(srcFile)[1]);
  }

  # Argument 'rows':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(rows)) {
    rows <- 1:nbrOfCells;
  } else {
    rows <- Arguments$getIndices(rows, range=c(1,nbrOfCells));
    rows <- sort(unique(rows));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing probe sequences");
  verbose && cat(verbose, "Pathname: ", getPathname(srcFile));
  verbose && cat(verbose, "Chip type: ", getChipType(srcFile));
  nbrOfColumns <- nbrOfColumns(this);
  verbose && cat(verbose, "Number of columns: ", nbrOfColumns);
  verbose && cat(verbose, "Rows:");
  verbose && str(verbose, rows);

  colClassPatterns <- c("probe(X|Y)Pos"="integer", "probeSequence"="character");

  count <- 0;
  CHUNK.SIZE <- as.integer(100e3);
  rowOffset <- as.integer(1);
  idxs <- 1:CHUNK.SIZE;
  while (length(rows) > 0) {
    if (length(rows) < CHUNK.SIZE) {
      idxs <- 1:length(rows);
    }

    rowsChunk <- rows[idxs];
    verbose && printf(verbose, "Row: %d-%d\n", min(rowsChunk), max(rowsChunk));
    df <- readDataFrame(srcFile, colClassPatterns=colClassPatterns, 
                        rows=rowsChunk, ..., verbose=less(verbose, 25));
    if (nrow(df) == 0)
      break;
    verbose && cat(verbose, "Data read:");
    verbose && str(verbose, df);

    if (nrow(df) < length(idxs)) {
      idxs <- 1:nrow(df);
    }

    cells <- nbrOfColumns*df[,"probeYPos"] + df[,"probeXPos"]+1;
    # Sanity check
    dups <- duplicated(cells);
    hasDuplicates <- any(dups);
    if (hasDuplicates) {
      setOfDups <- cells[dups];
      n <- length(setOfDups);
      verbose && print(verbose, head(df[dups,]));
      throw("Identified ", n, " duplicated cell indices: ", 
                                   paste(head(setOfDups), collapse=", "));
    }


    df[["probeXPos"]] <- df[["probeYPos"]] <- NULL;
    seqs <- df[,"probeSequence"];
    rm(df);
    gc <- gc();
    verbose && str(verbose, cells);
    verbose && str(verbose, seqs);

    updateSequenceStrings(this, cells=cells, seqs=seqs, 
                                               verbose=less(verbose, 25));
    count <- count + length(cells);

    # Next set of rows
    rows <- rows[-idxs];
  } # while(TRUE)

  verbose && exit(verbose);

  invisible(as.integer(count));
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getCdf", "AromaProbeSequenceTextFile", function(this, ...) {
  chipType <- getChipType(this, ...);
  AffymetrixCdfFile$byChipType(chipType);
})


setMethodS3("inferMmFromPmSequences", "AromaProbeSequenceTextFile", function(this, units=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


##   setMethodS3("dnaComplement", "character", function(s, ...) {
##     map <- c("A"="T", "C"="G", "G"="C", "T"="A");
##     s <- map[s];
##     s <- unname(s);
##     s;
##   }) # dnaComplement()
##   
##   setMethodS3("dnaComplement", "raw", function(s, ...) {
##     s <- rawToChar(s, multiple=TRUE);
##     dnaComplement(s);
##   }) # dnaComplement()
## 

  dnaComplement <- function(s, ...) {
    map <- c("A"="T", "C"="G", "G"="C", "T"="A");
    s <- map[s];
    s <- unname(s);
    s;
  } # dnaComplement()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  cdf <- getCdf(this);

  if (is.null(units)) {
    units <- seq(length=nbrOfUnits(cdf));
  } else {
    units <- Arguments$getIndices(units, range=c(1,nbrOfUnits(cdf)));
  }

  CHUNK.SIZE <- 100e3;
  nbrOfChunks <- ceiling(length(units) / CHUNK.SIZE);
  chunk <- 1;
  while (length(units) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    uu <- 1:min(length(units), CHUNK.SIZE);
    unitsChunk <- units[uu];

    # Read cell indices from CDF
    cells <- getCellIndices(cdf, units=unitsChunk, stratifyBy="pmmm", 
                                               unlist=TRUE, useNames=FALSE);

    # Create (PM,MM) pairs
    cells <- matrix(cells, nrow=2);
    rownames(cells) <- c("pm", "mm");
    verbose && cat(verbose, "(PM,MM) cell indices:");
    verbose && str(verbose, cells);

    # Read all PM sequences
    verbose && cat(verbose, "PM cell indices:");
    verbose && str(verbose, cells["pm",]);

    pmSeqs <- readSequenceMatrix(this, cells=cells["pm",], 
                                                verbose=less(verbose, 25));
    verbose && cat(verbose, "PM sequences:");
    verbose && str(verbose, pmSeqs);

    # Keep only (PM,MM) pairs for which we know the PM sequence
    keep <- which(!is.na(pmSeqs[,1]));
    verbose && printf(verbose, "Keeping %d of %d (%.1f%%) non-missing PM sequences\n", length(keep), nrow(pmSeqs), 100*length(keep)/nrow(pmSeqs));
    
    pmSeqs <- pmSeqs[keep,,drop=FALSE];
    cells <- cells["mm",keep];
    rm(keep);
    verbose && cat(verbose, "MM cell indices:");
    verbose && str(verbose, cells);

    # Complement 13th base for MM sequences
    mmSeqs <- pmSeqs;
    mmSeqs[,13] <- dnaComplement(mmSeqs[,13]);
    rm(pmSeqs);

    # Build MM sequences
    seqs <- seqMatrixToStrings(mmSeqs, ...);
    rm(mmSeqs);
    verbose && cat(verbose, "MM cell indices:");
    verbose && str(verbose, cells);
    verbose && cat(verbose, "MM sequences:");
    verbose && str(verbose, seqs);

#    o <- order(cells);
#    cells <- cells[o];
#    seqs <- seqs[o];
    updateSequenceStrings(this, cells=cells, seqs=seqs, 
                                                 verbose=less(verbose, 25));
    rm(cells, seqs);

    # Next chunk of units
    units <- units[-uu];
    chunk <- chunk + 1;
    verbose && exit(verbose);
  } # while()
})


setMethodS3("allocateFromCdf", "AromaProbeSequenceTextFile", function(static, ..., cdf) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  platform <- getPlatform(cdf);
  chipType <- getChipType(cdf);
  nbrOfRows <- nbrOfRows(cdf);
  nbrOfColumns <- nbrOfColumns(cdf);
  probeLengths <- 25;

  allocate(static, ..., platform=platform, chipType=chipType, nbrOfRows=nbrOfRows, nbrOfColumns=nbrOfColumns, probeLengths=probeLengths);
}, static=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2008-06-17
# o Added findByChipType() and byChipType().
# o Added allocateFromCdf().
# o Added countBases().
# o Major speed up of readSequenceStrings().
# 2008-06-16
# o Added fast isMissing().
# o Can now read and write probe sequences.
# 2008-06-14
# o Created.
############################################################################ 
