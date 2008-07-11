setMethodS3("readRawPairSequenceMatrix", "AromaProbeSequenceTextFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  seqs <- readRawSequenceMatrix(this, ..., map=as.raw(c(1:4, 0)), 
                                                       verbose=verbose);
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build paired sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pairSeqs <- matrix(as.raw(0), nrow=nrow(seqs), ncol=ncol(seqs)-1);

  # Only need to calculate non-missing sequences
  rr <- which(seqs[,1] != as.raw(0));
  if (length(rr) > 0) {
    for (cc in 1:ncol(pairSeqs)) {
      pair <- 4*as.integer(seqs[rr,cc]) + as.integer(seqs[rr,cc+1]);
      ok <- which(pair > 0);
      pair[ok] <- pair[ok] - as.integer(4);
      rm(ok);
      pair <- as.raw(pair);
      pairSeqs[rr,cc] <- pair;
      rm(pair);
    }
  }

  pairSeqs;
}, protected=TRUE)


setMethodS3("countPairs", "AromaProbeSequenceTextFile", function(this, cells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Read data
  pairSeqs <- readRawPairSequenceMatrix(this, ...);


  counts;
}, protected=TRUE)



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





############################################################################
# HISTORY:
# 2008-07-09
# o Added readRawPairSequenceMatrix().
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
