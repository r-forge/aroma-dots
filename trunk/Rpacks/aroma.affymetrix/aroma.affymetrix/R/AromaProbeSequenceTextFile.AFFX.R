# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
