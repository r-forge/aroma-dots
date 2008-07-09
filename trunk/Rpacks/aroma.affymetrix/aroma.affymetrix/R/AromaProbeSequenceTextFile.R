setConstructorS3("AromaProbeSequenceTextFile", function(...) {
  extend(AromaMicroarrayDataFile(...), "AromaProbeSequenceTextFile");
})


setMethodS3("allocate", "AromaProbeSequenceTextFile", function(this, name, tags=NULL, path=NULL, suffix=",probeSeqs.txt", platform, chipType, nbrOfRows, nbrOfColumns, probeLengths, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'name' and 'tags':
  fullname <- paste(c(name, tags), collapse=",");
  filename <- paste(fullname, suffix, sep="");
  pathname <- Arguments$getWritablePathname(filename, path=path, 
                                                  mustNotExist=!overwrite);

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
    nbrOfRows = nbrOfRows,
    nbrOfColumns = nbrOfColumns,
    probeLengths = probeLengths,
    probeSep = probeSep,
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


  verbose && cat(verbose, "Pathname: ", pathname);

  verbose && cat(verbose, "Cells to be interrogated:");
  verbose && str(verbose, cells);

  dataOffset <- getDataFileOffset(this);
  verbose && cat(verbose, "File offset to data section: ", dataOffset);

  seqLength <- as.integer(probeLengths+1);
  offset <- dataOffset + as.integer(1);
  if (is.null(cells)) {
    idxs <- seq(from=offset, to=offset+(seqLength*(nbrOfCells-1)), by=seqLength);
  } else {
    idxs <- offset + seqLength*cells;
  }
  rm(cells);

  verbose && cat(verbose, "File positions:");
  verbose && str(verbose, idxs);

  bfr <- readBinFragments(pathname, what="raw", idxs=idxs);
  rm(idxs);
  verbose && cat(verbose, "Raw data read:");
  verbose && str(verbose, bfr);

  isMissing <- (bfr == charToRaw(" "));
  rm(bfr);
  verbose && cat(verbose, "isMissing:");
  verbose && str(verbose, isMissing);
  verbose && str(verbose, which(isMissing));

  # Sanity check
  if (length(isMissing) != nbrOfCells) {
    throw("The length of the result vector does not match the number of cells interrogates: ", length(isMissing), " != ", nbrOfCells);
  }

  isMissing;
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
  if (is.null(cells)) {
    fromTo <- matrix(c(1,nbrOfCells), ncol=2);
  } else {
    fromTo <- seqToIntervals(cells);
  }

  # Translate to intervals in bytes
  seqLength <- as.integer(probeLengths+1);
  fromTo <- seqLength*(fromTo-1)+1;
  fromTo[,2] <- fromTo[,2] + probeLengths;
  fromTo <- dataOffset + fromTo;
  verbose && cat(verbose, "Number of byte intervals: ", nrow(fromTo));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reading raw data
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
  seqLength <- as.integer(probeLengths+1);
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
    unknown <- sprintf("%#02x (%s)", as.integer(unknown), 
                                        rawToChar(unknown, multiple=TRUE));
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

  # Memory clean up
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Ordering objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(cells)) {
    verbose && enter(verbose, "Ordering data");
    o <- order(cells);
    cells <- cells[o];
    gc <- gc();
    seqs <- matrix(seqs, nrow=seqLength);    
    seqs <- seqs[,o];
    seqs <- as.vector(seqs);
    rm(o);
    gc <- gc();
    verbose && print(verbose, gc);
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
    fromTo <- matrix(as.integer(c(1,nbrOfCells)), ncol=2);
  } else {
    fromTo <- seqToIntervals(cells);
  }
  verbose && str(verbose, fromTo);

  rm(cells);
  gc <- gc();
  verbose && print(verbose, gc);


  # Translate to intervals in bytes
  fromTo <- seqLength*(fromTo-as.integer(1))+as.integer(1);
  fromTo[,2] <- fromTo[,2] + probeLengths;
  fromTo <- dataOffset + fromTo;
  verbose && cat(verbose, "Number of byte intervals: ", nrow(fromTo));
  verbose && str(verbose, fromTo);

  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Writing raw data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  writeBinFragments(con=getPathname(this), seqs, idxs=fromTo);

  gc <- gc();
  verbose && print(verbose, gc);

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


setMethodS3("readSequences", "AromaProbeSequenceTextFile", function(this, ..., verbose=FALSE) {
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


setMethodS3("updateSequences", "AromaProbeSequenceTextFile", function(this, ..., seqs, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


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
  gc <- gc();

  seqs <- paste(seqs, collapse="\n");
  gc <- gc();
  seqs <- charToRaw(seqs);
  gc <- gc();
  seqs <- c(seqs, charToRaw("\n"));
  gc <- gc();
  verbose && print(verbose, gc);

  updateRawSequences(this, ..., seqs=seqs, verbose=less(verbose, 5));
})



setMethodS3("countBases", "AromaProbeSequenceTextFile", function(this, cells=NULL, ..., force=FALSE, verbose=FALSE) {
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



  verbose && enter(verbose, "Counting occurances of each base");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  key <- list();
  chipType <- getChipType(this);
  key <- list(method="countBases", class=class(this)[1], 
              chipType=chipType, fullname=getFullName(this), cells=cells);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    counts <- loadCache(key=key, dirs=dirs);
    if (!is.null(counts)) {
      verbose && cat(verbose, "Cached results found.");
      return(counts);
    }
  }
 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Optimize reading order?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(cells)) {
    cells <- 1:nbrOfCells;
    reorder <- FALSE;
  } else {
    o <- order(cells);
    cells <- cells[o];
    reorder <- TRUE;
    gc <- gc();
    verbose && print(verbose, gc);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Count in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  counts <- matrix(as.integer(NA), nrow=nbrOfCells, ncol=4);
  colnames(counts) <- c("A", "T", "C", "G");
  rawValues <- charToRaw(paste(colnames(counts), collapse=""));

  CHUNK.SIZE <- 1e6;
  offset <- 0;
  chunk <- 1;
  # Values to count/tabulate
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
    seqs <- readRawSequenceMatrix(this, cells=cellsChunk, ..., 
                                                verbose=less(verbose, 25));
    rm(cellsChunk);
    verbose && cat(verbose, "Sequences read:");
    verbose && str(verbose, seqs);

    # Read fewer cells?
    n2 <- nrow(seqs);
    if (n2 < n) {
      n <- n2;
      idxs <- 1:n;
    }

    # Count bases
    countsChunk <- matrixStats::rowTabulates(seqs, values=rawValues);

    # Identify missing probe sequences
    rowCounts <- rowSums(countsChunk);
    missing <- which(rowCounts == 0);
    rm(rowCounts);
    countsChunk[missing,] <- as.integer(NA);
    rm(missing);

    rm(seqs);
    verbose && cat(verbose, "Counts:");
    verbose && str(verbose, countsChunk);

    # Store
    counts[offset+idxs,] <- countsChunk;
    rm(countsChunk);
    gc <- gc();
    verbose && print(verbose, gc);

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

  # Cache results
  saveCache(key=key, dirs=dirs, counts);

  verbose && exit(verbose);

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


setMethodS3("importFromAffymetrixProbeTabFile", "AromaProbeSequenceTextFile", function(this, srcFile, rows=NULL, ..., ram=1, verbose=FALSE) {
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
                                            chipTypeSrc, " != ", chipType);
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

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(1e-3,Inf));

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
  CHUNK.SIZE <- as.integer(ram*1.5e6);
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

    cells <- nbrOfColumns*df[,"probeYPos"] + df[,"probeXPos"] + as.integer(1);

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
    rm(dups);

    # Clean up to save memory
    df[["probeXPos"]] <- df[["probeYPos"]] <- NULL;
    seqs <- df[,"probeSequence"];
    rm(df);
    gc <- gc();
    verbose && str(verbose, cells);
    verbose && str(verbose, seqs);

    updateSequences(this, cells=cells, seqs=seqs, verbose=less(verbose, 25));
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


setMethodS3("inferMmFromPmSequences", "AromaProbeSequenceTextFile", function(this, units=NULL, ..., safe=FALSE, ram=1, verbose=FALSE) {
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
  # Argument 'safe':
  safe <- Arguments$getLogical(safe);
  if (safe) {
    throw("Argument 'safe' was TRUE. There is no implementation available that infers MM cell indices safely from the CDF. The non-safe version assumes that stratifyBy=\"pmmm\" will return MMs in every 2nd index.");
  }

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(1e-3,Inf));

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

  CHUNK.SIZE <- as.integer(ram*100e3);
  nbrOfChunks <- ceiling(length(units) / CHUNK.SIZE);
  chunk <- 1;
  while (length(units) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    uu <- 1:min(length(units), CHUNK.SIZE);
    unitsChunk <- units[uu];
    verbose && cat(verbose, "Units in this chunk:");
    verbose && str(verbose, unitsChunk);

    # Read cell indices from CDF
    verbose && enter(verbose, "Querying CDF for cell indices");
    cells <- getCellIndices(cdf, units=unitsChunk, stratifyBy="pmmm", 
                                               unlist=TRUE, useNames=FALSE);
    rm(unitsChunk);
    # Sanity check
    if (length(cells) %% 2 != 0) {
      throw("Expected an even number of cell indices: ", length(cells));
    }
    verbose && exit(verbose);

    # Create (PM,MM) pairs
    cells <- matrix(cells, nrow=2);
    rownames(cells) <- c("pm", "mm");
    verbose && cat(verbose, "(PM,MM) cell indices:");
    verbose && str(verbose, cells);

    # Read all PM sequences
    verbose && cat(verbose, "PM cell indices:");
    verbose && str(verbose, cells["pm",]);

    pmSeqs <- readSequenceMatrix(this, cells=cells["pm",], 
                                                verbose=less(verbose, 5));
    verbose && cat(verbose, "PM sequences:");
    verbose && str(verbose, pmSeqs);
    gc <- gc();

    # Keep only (PM,MM) pairs for which we know the PM sequence
    keep <- pmSeqs[,1];
    keep <- is.na(keep);
    keep <- !keep;
    keep <- which(keep);
    verbose && printf(verbose, "Keeping %d of %d (%.1f%%) non-missing PM sequences\n", length(keep), nrow(pmSeqs), 100*length(keep)/nrow(pmSeqs));
    gc <- gc();
    verbose && print(verbose, gc);
    
    pmSeqs <- pmSeqs[keep,,drop=FALSE];
    gc <- gc();

    cells <- cells["mm",keep];
    rm(keep);
    verbose && cat(verbose, "MM cell indices:");
    verbose && str(verbose, cells);
    gc <- gc();

    # Complement 13th base for MM sequences
    mmSeqs <- pmSeqs;
    rm(pmSeqs);
    mmSeqs[,13] <- dnaComplement(mmSeqs[,13]);
    gc <- gc();
    verbose && print(verbose, gc);
    

    # Build MM sequences
    seqs <- seqMatrixToStrings(mmSeqs, ...);
    rm(mmSeqs);
    verbose && cat(verbose, "MM cell indices:");
    verbose && str(verbose, cells);
    verbose && cat(verbose, "MM sequences:");
    verbose && str(verbose, seqs);
    gc <- gc();
    verbose && print(verbose, gc);

#    o <- order(cells);
#    cells <- cells[o];
#    seqs <- seqs[o];
    updateSequences(this, cells=cells, seqs=seqs, verbose=less(verbose, 5));
    rm(cells, seqs);
    gc <- gc();
    verbose && print(verbose, gc);

    # Next chunk of units
    units <- units[-uu];
    rm(uu);

    chunk <- chunk + 1;
    verbose && exit(verbose);
  } # while()
})


setMethodS3("allocateFromCdf", "AromaProbeSequenceTextFile", function(static, cdf, ..., path=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Argument 'path':
  if (is.null(path)) {
    path <- filePath("annotationData", "chipTypes", 
                                         getChipType(cdf, fullname=FALSE));
  }


  platform <- getPlatform(cdf);
  chipType <- getChipType(cdf, fullname=FALSE);
  nbrOfRows <- nbrOfRows(cdf);
  nbrOfColumns <- nbrOfColumns(cdf);
  probeLengths <- 25;

  allocate(static, name=chipType, ..., path=path, platform=platform, chipType=chipType, nbrOfRows=nbrOfRows, nbrOfColumns=nbrOfColumns, probeLengths=probeLengths);
}, static=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("readSequenceStrings", "AromaProbeSequenceTextFile", function(this, ...) {
  readSequences(this, ...);
}, deprecated=TRUE)

setMethodS3("updateSequenceStrings", "AromaProbeSequenceTextFile", function(this, ...) {
  updateSequences(this, ...);
}, deprecated=TRUE)


############################################################################
# HISTORY:
# 2008-07-07
# o Now allocateFromCdf() only uses the 'name' part of the chip type.
# 2008-07-01
# o BUG FIX: isMissing(..., cells=c(1,n)) returned n values.
# o Now countBases() uses matrixStats::rowTabulates() to count bases instead
#   of matchprobes::basecounts(). The reason for this is that rowTabulates()
#   can work off the raw sequence matrices without having to translate them
#   to strings.  rowTabulates() is also not hardwired to count only ATCG:s,
#   which we will make use of later.
# 2008-06-21
# o Now countBases() caches results.
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
