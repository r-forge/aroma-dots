###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod computeAffinities
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
# }
#*/###########################################################################
setMethodS3("computeAffinities", "AffymetrixCdfFile", function(this, paths=NULL, force=FALSE, verbose=FALSE, ...) {
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

  chipType <- getChipType(this);


  verbose && enter(verbose, "Computing GCRMA probe affinities for ", nbrOfUnits(this), " units");

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


  # Checking cache
  key <- list(method="computeAffinities", class=class(this)[1], chipType=chipType);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res))
      return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate find probe sequence file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  verbose && enter(verbose, "Reading tab-delimited sequence file");
  sep <- "\t"

  # Note: Not all probe sequence tab files have column headers
  oLevel <- setDefaultLevel(verbose, -50);
  on.exit(setDefaultLevel(verbose, oLevel));

  verbose && enter(verbose, "Identifying first data row");
  lines <- readLines(psFile, n=2);
  hasHeader <- (regexpr("[Pp]robe.*[Ss]equence", lines[1]));
  if (hasHeader) {
    data <- lines[2];
  } else {
    data <- lines[1];
  }
  verbose && cat(verbose, "First data row: ", data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifing the X, Y, and sequence columns");
  data <- strsplit(data, split=sep)[[1]];
  nbrOfColumns <- length(data);

  # Get the unit name (assume first column)
  unitName <- data[1];

  # Find the sequence column(s)
  idxs <- grep(pattern="[actgACTG]{25}", data);
  if (length(idxs) == 0) {
    throw("Could not identify sequence column: ", paste(data, collapse=", "));
  } else if (length(idxs) > 1) {
    warning("Found more than one column containing sequence data. Using first."); 
  }
  seqcol <- idxs[1];

  # Identify potential (x,y) columns
  idxs <- grep(pattern="[0-9][0-9]*", data);
  if (length(idxs) < 2) {
    throw("Could not identify (x,y) columns: ", paste(data, collapse=", "));
  }
  values <- as.integer(data[idxs]);

  # Get the (x,y) CDF data for this unit
  unitNames <- getUnitNames(this);
  unitIdx <- match(unitName, unitNames);
  verbose && printf(verbose, "Unit: #%d (%s)\n", unitIdx, unitName);
  if (length(unitIdx) == 0)
    throw("Unit '", unitName, "' not found in CDF: ", getPathname(this));
  unitInfo <- readUnits(this, units=unitIdx);
  x <- applyCdfGroups(unitInfo, cdfGetFields, "x");
  x <- unlist(x, use.names=FALSE);
  y <- applyCdfGroups(unitInfo, cdfGetFields, "y");
  y <- unlist(y, use.names=FALSE);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  
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

  verbose && printf(verbose, "(X,Y,sequence) columns: (%d,%d,%d)\n", 
                                                      xcol, ycol, seqcol);
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  # read everybody in using read.table() with column classes specified
  # to speed things up
  colClasses <- rep("NULL", nbrOfColumns);
  colClasses[c(1,seqcol)] <- "character";
  colClasses[c(xcol,ycol)] <- "integer";
  names(colClasses)[c(1,xcol,ycol,seqcol)] <- c("name", "x", "y", "sequence");

  verbose && cat(verbose, "Column classes for read.table():");
  verbose && print(verbose, colClasses);
  sequenceInfo <- read.table(psFile, colClasses=colClasses, sep=sep, 
                                          skip=as.integer(hasHeader));
  colnames(sequenceInfo) <- names(colClasses)[colClasses != "NULL"];

  nbrOfSequences <- nrow(sequenceInfo);
  verbose && printf(verbose, "Loaded %d PM sequences", nbrOfSequences);
  verbose && str(verbose, sequenceInfo);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  setDefaultLevel(verbose, oLevel);

  # TODO: Reorder units according to CDF
  units <- match(sequenceInfo$name, unitNames);
  # Never used?!? /HB 2007-03-26
  rm(units, unitNames);
  verbose && exit(verbose);

  verbose && enter(verbose, "Get CDF cell indices for all units");
  indices <- getCellIndices(this, useNames=FALSE, unlist=TRUE);
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
 
# assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  indexPm <- sequenceInfo$y * dimension[1] + sequenceInfo$x + 1
  indexMmPutative <- indexPm + dimension[1];


# match up putative MM with actual list of MMs from CDF
  matches <- match(indexMmPutative, indices[isMm]);
  indexMm <- indices[isMm][matches];
  notNA <- which(!is.na(indexMm));
  indexMm <- indexMm[notNA];

  rm(indexMmPutative, matches); # Not needed anymore

# now calculate affinities - code reused from compute.affinities() in
# gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  # To please R CMD check on R v2.6.0
  affinity.spline.coefs <- NULL; rm(affinity.spline.coefs);
  data("affinity.spline.coefs"); # A tiny object from 'gcrma'.

  affinity.basis.matrix <- splines::ns(1:25, df=length(affinity.spline.coefs)/3);

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5]);
  T13 <- 0;
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10]);
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15]);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  apm <- vector("double", nbrOfSequences);
  amm <- vector("double", nbrOfSequences);

  if (verbose) {
    cat(verbose, "Progress (counting to 100): ");
    pb <- ProgressBar(stepLength=100/(nbrOfSequences/1000));
    reset(pb);
  }

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

  # create a vector to hold affinities and assign values to the appropriate
  # location in the vector
  affinities <- rep(NA, dimension[1]*dimension[2]);
  affinities[indexPm] <- apm;
  affinities[indexMm] <- amm[notNA];
  verbose && exit(verbose);
  rm(dimension, indexPm, indexMm, apm, amm, notNA); # Not needed anymore

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
