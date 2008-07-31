###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod computeAffinities2
#
# @title "Calculates probe affinities from sequence"
#
# \description{
#  @get "title".
#
# Adapted from @see "gcrma::compute.affinities" in the \pkg{gcrma} package.
# Attempts to find the tab-separated probe sequence file associated with
# a particular CDF, and matches sequence to probe index in order to assign
# an affinity to each probe.
# }
#
# @synopsis
#
# \arguments{
#   \item{paths}{A @character variable containing location(s) to look for
#    the probe sequence file.  Multiple locations should be separated by
#    semi-colons.}
#   \item{force}{If @FALSE, cached results is returned, if available.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{safe}{If @TRUE and the probe-tab file has no column names, then
#    extra care is take to make sure the correct data columns are extracted.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector of (log2) probe affinities, of length equal
#  to the total number of features on the array.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
#   Modified by Henrik Bengtsson.
# }
#*/###########################################################################
setMethodS3("computeAffinities2", "AffymetrixCdfFile", function(this, paths=NULL, safe=TRUE, force=FALSE, verbose=FALSE, ...) {
  # Try to load all required package first
  require("gcrma", quietly=TRUE) || throw("Package not loaded: gcrma");
  require("matchprobes", quietly=TRUE) || throw("Package not loaded: matchprobes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Computing GCRMA probe affinities for ", nbrOfUnits(this), " units");

  # Get the chip type
  chipTypeFull <- getChipType(this, fullname=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Looking for MMs (and PMs) in the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identify PMs and MMs among the CDF cell indices");
  isPm <- isPm(this);
  verbose && str(verbose, isPm);
  verbose && summary(verbose, isPm);
  nbrOfPMs <- sum(isPm);
  verbose && cat(verbose, "MMs are defined as non-PMs");
  isMm <- !isPm;
  nbrOfMMs <- sum(isMm);
  verbose && cat(verbose, "Number of PMs: ", nbrOfPMs);
  verbose && cat(verbose, "Number of MMs: ", nbrOfMMs);
  verbose && exit(verbose);

  # Sanity check
  if (nbrOfMMs == 0) {
#    throw("Cannot calculate gcRMA probe affinities. The CDF contains no MMs: ", chipTypeFull);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking for cache results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="computeAffinities", class=class(this)[1], chipTypeFull=chipTypeFull);
  dirs <- c("aroma.affymetrix", chipTypeFull);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res))
      return(res);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate find probe sequence file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The probe-sequence does not depend on the CDF but only the chip type,
  # which is why we ignore any tags for the CDF.
  chipType <- getChipType(this, fullname=FALSE);
  psFile <- AffymetrixProbeTabFile$findByChipType(chipType=chipType, 
                                         paths=paths, verbose=less(verbose));
  if (is.null(psFile))
    throw("Could not locate probe sequence file for chip type: ", chipType);

  verbose && cat(verbose, "Pathname of located probe-tab file: ", psFile);

  
# read in probe sequence data

  dimension <- getDimension(this);

#  con <- file(psFile, open="r")
#  on.exit({
#    if (!is.null(con))
#      close(con);
#  });


  # Things needed below
  unitNames <- getUnitNames(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying columns for (x,y, sequence)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying columns for (x,y, sequence)");

  verbose && enter(verbose, "Reading tab-delimited sequence file");
  sep <- "\t";
  # Note: Not all probe sequence tab files have column headers
  oLevel <- setDefaultLevel(verbose, -10);
  on.exit(setDefaultLevel(verbose, oLevel));
  lines <- readLines(psFile, n=2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Checking if the file contains column names or not");
  hasHeader <- (regexpr("[Pp]robe.*[Ss]equence", lines[1]) != -1);
  verbose && cat(verbose, "Probe-tab file has header: ", hasHeader);
  verbose && exit(verbose);

#  hasHeader <- FALSE; lines <- lines[-1]; # Fake it for now.
  if (hasHeader) {
    verbose && enter(verbose, "Identifying columns by column names");
    columnNames <- strsplit(lines[1], split=sep)[[1]];
    nbrOfColumns <- length(columnNames);
    verbose && cat(verbose, "Column names: ", paste(columnNames, collapse=", "));
    columnNames <- toCamelCase(columnNames);
    verbose && cat(verbose, "Camel-case column names: ", paste(columnNames, collapse=", "));

    unitColumn <- grep("probeSetID", columnNames);
    seqcol <- grep("probe.*[Ss]equence", columnNames);
    xcol <- grep("probeX", columnNames);
    ycol <- grep("probeY", columnNames);
    data <- lines[2];    
    data <- strsplit(data, split=sep)[[1]];
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Identifying columns by guessing from data");
    data <- lines[1];
    verbose && cat(verbose, "First data row: ", data);
    data <- strsplit(data, split=sep)[[1]];
    nbrOfColumns <- length(data);
    verbose && cat(verbose, "Number of column: ", nbrOfColumns);

    verbose && enter(verbose, "Identifying the column contain unit names");
    for (unitColumn in 1:2) {
      unitName <- data[unitColumn];
      unitIdx <- match(unitName, unitNames);
      if (!is.na(unitIdx))
        break;
    }

#    unitIdx <- 1;  # Fake it for now.

    if (is.na(unitIdx)) {
      throw("Failed to locate column containing unit names. Neither of column one and two contain a valid unit name according to the CDF.");
    }

    verbose && cat(verbose, "Column: ", unitColumn);
    verbose && cat(verbose, "Example (unit, unitName): (", unitIdx, ", '", unitName, "')");
    verbose && exit(verbose);

    verbose && enter(verbose, "Identifying the column containing sequences");
    # Find the sequence column(s)
    idxs <- grep(pattern="[actgACTG]{25}", data);
    if (length(idxs) == 0) {
      throw("Could not identify sequence column: ", paste(data, collapse=", "));
    } else if (length(idxs) > 1) {
      msg <- paste("Found more than one column containing sequence data. Using first: ", paste(idxs, collapse=", "), sep=""); 
      verbose && cat(verbose, msg);
      warning(msg);
    }
    seqcol <- idxs[1];
    verbose && cat(verbose, "Column: ", seqcol);
    verbose && cat(verbose, "Example: ", data[seqcol]);
    verbose && exit(verbose);

    verbose && enter(verbose, "Identifying the column containing (x,y)");
    # Identify potential (x,y) columns
    idxs <- grep(pattern="[0-9][0-9]*", data);
    
    if (length(idxs) < 2) {
      throw("Could not identify (x,y) columns: ", paste(data, collapse=", "));
    } else if (length(idxs) > 2) {
      # We assume that it starts after the unit name, so exclude anything 
      # before that.
      idxs <- idxs[idxs > unitColumn];
  
      if (length(idxs) > 2) {
        msg <- paste("Found more than two columns containing possible (x,y) data. Using the first two: ", paste(idxs, collapse=", "), sep=""); 
        verbose && cat(verbose, msg);
        warning(msg);
        idxs <- idxs[1:2];
      }
    }
    verbose && cat(verbose, "Unordered (x,y) columns: ", paste(idxs, collapse=", "));
    verbose && cat(verbose, "Example: (", paste(data[idxs], collapse=","), ")");
    verbose && exit(verbose);

    if (safe) {
      verbose && enter(verbose, "Identifying the order of the (x,y) columns");
      values <- as.integer(data[idxs]);
    
      verbose && enter(verbose, "Reading all (x,y) from the CDF for the unit");
      verbose && cat(verbose, "Unit: ", unitIdx);
      # Get the (x,y) CDF data for this unit
      unitInfo <- readUnits(this, units=unitIdx);
      x <- applyCdfGroups(unitInfo, cdfGetFields, "x");
      x <- unlist(x, use.names=FALSE);
      y <- applyCdfGroups(unitInfo, cdfGetFields, "y");
      y <- unlist(y, use.names=FALSE);
      verbose && str(verbose, x);
      verbose && str(verbose, y);
      verbose && exit(verbose);

      # Scan for x coordinate
      idxs <- which(values %in% x);
      if (length(idxs) == 0)
        throw("Could not identify X column: ", paste(data, collapse=", "));
      # If more than one match, assume first
      xcol <- idxs[1];
      ycol <- xcol + 1;
      # Get the (x,y) coordinate (in the sequence file)
      xySeq <- c(values[xcol], values[ycol]);
      verbose && cat(verbose, "(X,Y) in sequence file:");
      verbose && print(verbose, xySeq);

      # Now, find the same in the CDF file
      kk <- which(xySeq[1] == x);
      xyCdf <- matrix(c(x[kk],y[kk]), ncol=2);
      verbose && cat(verbose, "(X,Y):s in CDF file with matching X coordinate:");
      verbose && print(verbose, xyCdf);
      rr <- which(apply(xyCdf, MARGIN=1, FUN=identical, xySeq));
      if (length(rr) == 0)
        throw("Could not identify matching (X,Y) in CDF.");
      xyCdf <- xyCdf[rr,];
    
      rm(xySeq, xyCdf, rr); # Never used?!? /HB 2007-03-26
  
      verbose && exit(verbose);
    } else {
      xcol <- idxs[1];
      ycol <- idxs[2];
    }

    verbose && cat(verbose, "Final (x,y) columns: ", xcol, ", ", ycol);
    verbose && cat(verbose, "Example: (", data[xcol], ",", data[ycol], ")");

    verbose && exit(verbose);
  } # if (hasHeader)

  idxs <- c(xcol, ycol, seqcol);
  verbose && cat(verbose, "Identified (x,y, sequence) columns: ", paste(idxs, collapse=", "));
  verbose && cat(verbose, "Example: (x,y, sequence)=(", paste(data[idxs], collapse=", "), ")");
  if (length(idxs) != 3)
    throw("Failed to locate all columns");


  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading (x,y,sequence) data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading (x, y, sequence) data");
  # read everybody in using read.table() with column classes specified
  # to speed things up
  verbose && enter(verbose, "Setting up colClasses");
  colClasses <- rep("NULL", nbrOfColumns);
  colClasses[c(xcol,ycol)] <- "integer";
  colClasses[seqcol] <- "character";
  names(colClasses)[c(xcol,ycol,seqcol)] <- c("x", "y", "sequence");
  verbose && cat(verbose, "Column classes for read.table():");
  verbose && print(verbose, colClasses);
  verbose && exit(verbose);

  skip <- as.integer(hasHeader);
  verbose && cat(verbose, "skip: ", skip);
  verbose && cat(verbose, "Additional arguments:");
  verbose && str(verbose, list(...));
  data <- read.table(psFile, colClasses=colClasses, sep=sep, skip=skip, ...);
  colnames(data) <- names(colClasses)[colClasses != "NULL"];

  nbrOfSequences <- nrow(data);
  verbose && printf(verbose, "Loaded %d PM sequences", nbrOfSequences);
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  setDefaultLevel(verbose, oLevel);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying sequences for which there exists a (PM,MM) pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying sequences for which there exists a (PM,MM) pairs");

  verbose && enter(verbose, "Calculating the cell indices from (x,y)");
  cells <- data$y * dimension[1] + data$x + as.integer(1);
  cells <- as.integer(cells);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  verbose && exit(verbose);

  verbose && enter(verbose, "Get sequences");
  seqs <- data$sequence;
  verbose && cat(verbose, "Sequences:");
  verbose && str(verbose, seqs);
  verbose && exit(verbose);

  # Not needed anymore
  rm(data);
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


  verbose && enter(verbose, "Getting cell indices for all units in the CDF");
  cdfCells <- getCellIndices(this, useNames=FALSE, unlist=TRUE);
  verbose && str(verbose, cdfCells);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


# probe sequence tab-separated files generally only contain the PM
# sequence (since we know that MM sequence always has base 13 complemented).
# This is true for expression arrays, but possibly not for genotyping arrays.
# We will calculate affinity for PM and MM anyway, and then
# try to match MM indices from the CDF file with their appropriate PMs (and
# hence assign the correct affinity).

  verbose && enter(verbose, "Identifying cell indices for MMs");

  verbose && cat(verbose, "MM indices in the CDF:");
  cdfCellsMM <- cells[isMm];
  verbose && str(verbose, cdfCellsMM);
  rm(cdfCells);

  # Assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  cellsMmPutative <- cells + dimension[1];
  verbose && cat(verbose, "MM indices in the probe-tab file (assuming one row below PMs):");
  verbose && str(verbose, cellsMmPutative);

  verbose && cat(verbose, "MM indices in the probe-tab file (after matching to CDF):");

  
  verbose && enter(verbose, "Filtering out PM sequences for which there is no MM probe in the CDF");
  # match up putative MM with actual list of MMs from CDF
  pmHasMm <- (cellsMmPutative %in% cdfCellsMM);
  rm(cellsMmPutative, cdfCellsMM); # Not needed anymore
  verbose && summary(verbose, pmHasMm);

  cells <- cells[pmHasMm];
  seqs <- seqs[pmHasMm];
  verbose && summary(verbose, pmHasMm);
  rm(pmHasMm);
  verbose && exit(verbose);

  verbose && cat(verbose, "Indices of PM cells for existing (PM,MM) pairs with sequence information");
  verbose && str(verbose, cells);

  nbrOfSequences <- length(cells);

  # Sanity check
  if (nbrOfSequences == 0) {
#    throw("Cannot calculate gcRMA probe affinities. The probe-tab contains no sequences of PMs which a matching MM: ", chipTypeFull);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate probe-sequence affinities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # now calculate affinities - code reused from compute.affinities() in
  # gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  verbose && enter(verbose, "Setting up prior probe-sequence parameters");

  # To please R CMD check on R v2.6.0
  affinity.spline.coefs <- NULL; rm(affinity.spline.coefs);
  data("affinity.spline.coefs"); # A tiny object from 'gcrma'.

  affinity.basis.matrix <- splines::ns(1:25, df=length(affinity.spline.coefs)/3);
  verbose && cat(verbose, "affinity.basis.matrix:");
  verbose && str(verbose, affinity.basis.matrix);

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5]);
  T13 <- 0;
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10]);
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15]);

  verbose && cat(verbose, "(A13,T13,C13,G13): ", paste(c(A13,T13,C13,G13), collapse=", "));

  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Allocating probe-affinity vectors for (PM,MM) pairs");
  verbose && cat(verbose, "Number of sequences: ", nbrOfSequences);
  apm <- vector("double", nbrOfSequences);
  amm <- vector("double", nbrOfSequences);
  verbose && exit(verbose);


  verbose && enter(verbose, "Calculating affinities");

  if (verbose) {
    cat(verbose, "Progress (counting to 100): ");
    pb <- ProgressBar(stepLength=100/(nbrOfSequences/1000));
    reset(pb);
  }

  # Speed up
  T13A13 <- T13 - A13;
  C13G13 <- C13 - G13;

  for (ii in seq(along = apm)) {
    if (verbose && ii %% 1000 == 0)
      increase(pb);
    # Get a 4x25 matrix with rows A, C, G, and T.
    charMtrx <- .Call("gcrma_getSeq", sequenceInfo$sequence[ii], 
                                                       PACKAGE="gcrma");

    A <- cbind(charMtrx[1, ] %*% affinity.basis.matrix, 
               charMtrx[2, ] %*% affinity.basis.matrix, 
               charMtrx[3, ] %*% affinity.basis.matrix);
    apm[ii] <- A %*% affinity.spline.coefs;
    if (charMtrx[1, 13] == 1) {
      amm[ii] <- apm[ii] + T13A13;  # + T13 - A13
    } else {
      if (charMtrx[4, 13] == 1) {
        amm[ii] <- apm[ii] - T13A13;  # + A13 - T13
      } else {
        if (charMtrx[3, 13]) {
          amm[ii] <- apm[ii] + C13G13;  # + C13 - G13
        } else {
          amm[ii] <- apm[ii] - C13G13;  # + G13 - C13
        }
      }
    }
  } # for (ii in ...)
  rm(charMtrx, A); # Not needed anymore
  verbose && cat(verbose, "");
  verbose && exit(verbose);


  # create a vector to hold affinities and assign values to the appropriate
  # location in the vector
  naValue <- as.double(NA);
  affinities <- rep(naValue, dimension[1]*dimension[2]);
  affinities[cells] <- apm;
  affinities[cells+dimension[1]] <- amm;
  rm(apm, amm, cells);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Saving to cache
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, affinities, comment=comment, dirs=dirs);

  verbose && exit(verbose);
  
  affinities;
}, private=TRUE)





############################################################################
# HISTORY:
# 2008-07-30
# o UPDATE: Now computeAffinities() for AffymetrixCdfFile only computes
#   affinities for (PM,MM) pairs (according to the CDF) for which there
#   exists a probe sequence for the PM in the probe-tab file.
# o BUG FIX: computeAffinities() for AffymetrixCdfFile searched for the
#   probe-tab file using the chip type given by the fullname of the CDF
#   and not the basic name.
# 2007-09-06
# o Made computeAffinities() more memory efficient since it is using the
#   new getCellIndices(cdf, useNames=FALSE, unlist=TRUE).
# 2007-06-07
# o BUG FIX: If an Affymetrix probe tab file is not found for the chip type,
#   computeAffinitities() of AffymetrixCdfFile would throw "Error in 
#   paste(..., sep = sep, collapse = collapse): object "pattern" not found"
#   instead of an intended and more information error.
# 2007-02-06
# o Now computeAffinities() is locating the Affymetrix' probe-sequence file
#   using AffymetrixProbeTabFile$findByChipType().
# 2006-10-09
# o Added caching to file.
# o Now default probe affinity is NA (for non estimated probes).
# o Added a progress bar to the calculations.
# o Made internal reader smart so it guess the X, Y, and sequence columns
#   by reading the data and comparing to CDF information.
# o Added verbose statements.
# o Made sure connection is closed regardless how the method exits.
# 2006-10-04
# o Debugged, tested, docs/comments added.
# 2006-10-01
# o Created.
############################################################################
