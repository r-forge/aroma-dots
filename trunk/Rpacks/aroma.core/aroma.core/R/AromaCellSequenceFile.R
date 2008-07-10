setConstructorS3("AromaCellSequenceFile", function(...) {
  extend(AromaCellTabularBinaryFile(...), "AromaCellSequenceFile");
})

setMethodS3("getFilenameExtension", "AromaCellSequenceFile", function(static, ...) {
  "acs";
}, static=TRUE)


setMethodS3("getColumnNames", "AromaCellSequenceFile", function(this, ...) {
  c(sprintf("b%02d", seq(from=1, to=nbrOfColumns(this)-1)), "targetStrand");
})


setMethodS3("getProbeLength", "AromaCellSequenceFile", function(this, ...) {
  as.integer(nbrOfColumns(this, ...) - 1);
})


setMethodS3("readSequenceMatrix", "AromaCellSequenceFile", function(this, cells=NULL, positions=seq(length=getProbeLength(this)), drop=FALSE, what=c("raw", "character"), naValue=NA, ...) {
  # Argument 'positions':
  positions <- Arguments$getIndices(positions, range=c(1, getProbeLength(this)));

  # Argument 'what':
  what <- match.arg(what);

  # Read data
  res <- this[cells, positions];

  # The raw to character map  
  map <- as.raw(0:4);
  names(map) <- c(" ", "A", "C", "G", "T");

  # Flatten (data frame)
  dim <- dim(res);
  res <- unlist(res, use.names=FALSE);

  # Coerce to character strings?
  if (what == "character") {
    res <- as.integer(res) + as.integer(1);
    names <- names(map);
    names[names == " "] <- naValue;
    names(map) <- names;
    res <- names[res];
  }

  # Coerce to matrix
  dim(res) <- dim;

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  # Add 'map' attribute
  attr(res, "map") <- map;

  res;
})

setMethodS3("readPairSequenceMatrix", "AromaCellSequenceFile", function(this, ...) {
  readNeighborSequenceMatrix(this, nbrOfNeighbors=2, ...);
})


setMethodS3("readNeighborSequenceMatrix", "AromaCellSequenceFile", function(this, ..., nbrOfNeighbors, drop=FALSE, what=c("character", "raw"), naValue=NA, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'nbrOfNeighbors':
  nbrOfNeighbors <- Arguments$getInteger(nbrOfNeighbors, range=c(1, getProbeLength(this)));

  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading sequence matrix");
  seqWhat <- "raw";
  if (nbrOfNeighbors == 1)
    seqWhat <- what;
  seqs <- readSequenceMatrix(this, ..., what=seqWhat, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Nothing more to do?
  if (nbrOfNeighbors == 1) {
    return(seqs);
  }

  # Identify nucleotides
  map <- attr(seqs, "map");
  names <- names(map);
  names <- names[names != " "];

  neighborNames <- names;
  for (kk in seq(from=2, to=nbrOfNeighbors)) {
    neighborNames <- outer(neighborNames, names, FUN=paste, sep="");
  }
  neighborNames <- as.vector(neighborNames);
  neighborNames <- sort(neighborNames);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build neighbored sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating return matrix");
  res <- matrix(as.raw(0), nrow=nrow(seqs), ncol=ncol(seqs)-1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying non-missing sequences");
  # Only need to calculate non-missing sequences
  rr <- which(seqs[,1] != as.raw(0));
  verbose && exit(verbose);

  verbose && enter(verbose, "Mapping sequences to neighbor-group sequences");
  basis <- rev(as.integer(4^(1:nbrOfNeighbors-1)));
  if (length(rr) > 0) {
    for (cc in 1:ncol(res)) {
      values <- as.integer(seqs[rr,cc]);
      values <- values - as.integer(1);
      values <- basis[1] * values;
      neighbors <- values;
      rm(values);

      for (tt in 2:nbrOfNeighbors) {
        values <- as.integer(seqs[rr,cc+tt-1]);
        if (tt < nbrOfNeighbors) {
          values <- values - as.integer(1);
          values <- basis[tt] * values;
        }
        neighbors <- neighbors + values;
        rm(values);
      }

      neighbors <- as.raw(neighbors);
      res[rr,cc] <- neighbors;
      rm(neighbors);
    } # for (cc ...)
  } 
  verbose && exit(verbose);

  map <- as.raw(0:length(neighborNames));
  names(map) <- c(" ", neighborNames);

  # Coerce to character strings?
  if (what == "character") {
    dim <- dim(res);
    res <- as.integer(res) + as.integer(1);
    names <- names(map);
    names[names == " "] <- naValue;
    names(map) <- names;
    res <- names[res];
    dim(res) <- dim;
  }

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  attr(res, "map") <- map;

  res;
}, protected=TRUE)



setMethodS3("readSequences", "AromaCellSequenceFile", function(this, ..., naValue=NA) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'naValue':


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- readSequenceMatrix(this, ..., what="raw");
  map <- attr(res, "map");

  nbrOfCells <- nrow(res);
  nbrOfPositions <- ncol(res);

  # Setup naValue
  if (is.na(naValue)) {
  } else {
    naValue <- as.character(naValue);
    naValue <- rep(naValue, nbrOfPositions);
    naValue <- paste(naValue, collapse="");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup and allocation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify non-missing sequences
  idxs <- (res[,1] != map[" "]);
  idxs <- which(idxs);

  # Keep only those
  res <- res[idxs,,drop=FALSE];

  # Allocate return vector with missing values set
  seqs <- rep(naValue, times=nbrOfCells);
  seqs[idxs] <- "";  # Redo non-missing

  # Nothing more to do?
  if (nrow(res) == 0) {
    return(seqs);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remap
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  values <- names(map);
  values <- paste(values, collapse="");
  values <- charToRaw(values);

  for (kk in seq(along=map)) {
    idxsT <- (res == map[kk]);
    idxsT <- which(idxsT);
    res[idxsT] <- values[kk];
  }
  rm(idxsT);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build sequence strings
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while (ncol(res) > 0) {
    bases <- rawToChar(res[,1], multiple=TRUE);
    seqs[idxs] <- paste(seqs[idxs], bases, sep="");
    res <- res[,-1,drop=FALSE];
  }
  gc <- gc();

  seqs;
})


setMethodS3("isMissing", "AromaCellSequenceFile", function(this, ...) {
  res <- readSequenceMatrix(this, ..., positions=1, what="raw", drop=TRUE);
  res <- (res == as.raw(0));
  res;
}, protected=TRUE)


setMethodS3("countBases", "AromaCellSequenceFile", function(this, bases=c("A", "C", "G", "T"), drop=FALSE, ...) {
  # Tabular nucleotides
  counts <- countBasesInternal(this, ...);

  # Identify missing sequences
  isMissing <- which(counts[,1] > 0);

  # Keep only bases of interest
  counts <- counts[,bases,drop=FALSE];

  # Set missing values
  counts[isMissing,] <- NA;

  # Drop singleton dimensions?
  if (drop) {
    counts <- drop(counts);
  }

  counts;
})

setMethodS3("countBasesInternal", "AromaCellSequenceFile", function(this, cells=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
    nbrOfCells <- length(cells);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Counting occurances of each nucleotide");

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
  # Sum over positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  counts <- matrix(0, nrow=nbrOfCells, ncol=5);

  map <- NULL;
  nbrOfPositions <- getProbeLength(this);
  for (pp in seq(length=nbrOfPositions)) {
    verbose && enter(verbose, sprintf("Position #%d of %d", pp, nbrOfPositions));
    seqs <- readSequenceMatrix(this, cells=cells, positions=pp, what="raw", drop=TRUE);

    if (is.null(map)) {
      map <- attr(seqs, "map");
      colnames(counts) <- names(map);
    }

    # Add to counts
    for (bb in seq(length=ncol(counts))) {
      idxs <- (seqs == map[bb]);
      idxs <- which(idxs);
      counts[idxs,bb] <- counts[idxs,bb] + as.integer(1);
      rm(idxs);
    }
    rm(seqs);
    
    verbose && exit(verbose);
  } # for (pp ...)

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
 

setMethodS3("allocate", "AromaCellSequenceFile", function(static, ..., platform, chipType, footer=list()) {
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'footer':
  if (is.null(footer)) {
  } else if (!is.list(footer)) {
    throw("Argument 'footer' must be NULL or a list: ", class(footer)[1]);
  }

  footer <- c(
    list(
      createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
      platform=platform,
      chipType=chipType
    ), 
    footer
  );

  probeLengths <- 25;
  nbrOfColumns <- probeLengths+1;
str(111);
  res <- allocate.AromaMicroarrayTabularBinaryFile(static, ..., 
                 types=rep("raw", nbrOfColumns), sizes=rep(1,nbrOfColumns), 
                                                            footer=footer);
str(222);

  res;
}, static=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethodS3("allocateFromCdf", "AromaCellSequenceFile", function(static, cdf, path=getPath(cdf), tags=NULL, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(cdf);

  # Exclude 'monocell' tags (AD HOC)
  chipType <- gsub(",monocell", "", chipType);

  # Get platform
  platform <- getPlatform(cdf);

  # Number of cells
  nbrOfCells <- nbrOfCells(cdf);

  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create microarray tabular binary file
  allocate(static, filename=filename, path=path, nbrOfRows=nbrOfCells, 
                                platform=platform, chipType=chipType, ...);
}, static=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethodS3("importFromAromaProbeSequenceTextFile", "AromaCellSequenceFile", function(this, srcFile, cells=seq(length=nbrOfCells(this)), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'srcFile':
  if (!inherits(srcFile, "AromaProbeSequenceTextFile")) {
    throw("Argument 'srcFile' is not an AromaProbeSequenceTextFile: ", 
                                                             class(srcFile)[1]);
  }

  # Argument 'cells':
  cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells(this)));
  cells <- unique(cells);
  cells <- sort(cells);


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert compatibility
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (getPlatform(srcFile) != getPlatform(this)) {
    throw("The platform of argument 'srcFile' is not '", getPlatform(this), "': ", getPlatform(srcFile));
  }

  if (getChipType(srcFile, fullname=FALSE) != getChipType(this, fullname=FALSE)) {
    throw("The chip type of argument 'srcFile' is not '", getChipType(this), "': ", getChipType(srcFile));
  }

  if (nbrOfCells(srcFile) != nbrOfCells(this)) {
    throw("The number of cells of argument 'srcFile' is not ", nbrOfCells(this), ": ", nbrOfCells(srcFile));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Import data in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  positions <- 1:getProbeLength(this);
  CHUNK.SIZE <- 500e3;
  cc <- 1:CHUNK.SIZE;
  while(length(cells) > 0) {
    verbose && cat(verbose, "Cells left to import:");
    verbose && str(verbose, cells);
    if (length(cells) < CHUNK.SIZE) {
      cc <- 1:length(cells);
    }
    cellsChunk <- cells[cc];
    verbose && cat(verbose, "Cells to read in this chunk:");
    verbose && str(verbose, cellsChunk);

    # Read data
    seqs <- readRawSequenceMatrix(srcFile, cells=cellsChunk, 
                           map=as.raw(c(1:4,0)), verbose=less(verbose, 5));
    verbose && cat(verbose, "Raw sequence matrix read:");
    verbose && str(verbose, seqs);

    # Sanity check
    if (nrow(seqs) != length(cellsChunk)) {
      throw("Internal error.");
    }
    if (ncol(seqs) != length(positions)) {
      throw("Internal error.");
    }

    # Write data
    this[cellsChunk,positions] <- seqs;
    rm(seqs, cellsChunk);

    # Next chunk
    cells <- cells[-cc];

    gc <- gc();
  } # while(...)

  invisible(this);
}, protected=TRUE)

############################################################################
# HISTORY:
# 2008-07-09
# o Created from AromaUgpFile.R.
############################################################################
