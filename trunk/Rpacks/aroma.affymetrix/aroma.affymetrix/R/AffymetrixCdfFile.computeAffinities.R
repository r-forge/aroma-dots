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
  require(gcrma, quietly=TRUE) || throw("Package not loaded: gcrma");
  require(splines, quietly=TRUE) || throw("Package not loaded: splines");
  require(matchprobes, quietly=TRUE) || throw("Package not loaded: matchprobes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  chipType <- getChipType(this)


  verbose && enter(verbose, "Computing GCRMA probe affinities for ", nbrOfUnits(this), " units");

  # Check cache
  key <- list(method="computeAffinities", class=class(this)[1], chipType=chipType);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res))
      return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # try to find probe sequence file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for probe sequence file");
  if (is.null(paths)) {
    paths <- paste(".",
                   getOption("AFFX_SEQUENCE_PATH"),
                   Sys.getenv("AFFX_SEQUENCE_PATH"),
                   "sequence/", "data/sequence/",
                   getOption("AFFX_CDF_PATH"),
                   Sys.getenv("AFFX_CDF_PATH"),
                   "cdf/", "data/cdf/",
                   sep=";", collapse=";");
  }

  pattern <- paste(chipType, "_probe_tab", sep="");
  verbose && cat(verbose, "Filename pattern: ", pattern);

  psFile <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
  if (is.null(psFile)) 
    throw("Could not locate probe sequence file: ", pattern);
  verbose && cat(verbose, "Pathname: ", psFile);
  verbose && exit(verbose);

  
# read in probe sequence data

  dimension <- getDimension(this)

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
  y <- applyCdfGroups(unitInfo, cdfGetFields, "y");
  x <- unlist(x, use.names=FALSE);
  y <- unlist(y, use.names=FALSE);
  
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

  verbose && printf(verbose, "(X,Y,sequence) columns: (%d,%d,%d)\n", 
                                                      xcol, ycol, seqcol);
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

  setDefaultLevel(verbose, oLevel);

  # TODO: Reorder units according to CDF
  units <- match(sequenceInfo$name, unitNames);
  verbose && exit(verbose);

  verbose && enter(verbose, "Get CDF cell indices for all units");
  indices <- unlist(getCellIndices(this), use.names=FALSE);
  verbose && exit(verbose);


# probe sequence tab-separated files generally only contain the PM
# sequence (since we know that MM sequence always has base 13 complemented).
# This is true for expression arrays, but possibly not for genotyping arrays.
# We will calculate affinity for PM and MM anyway, and then
# try to match MM indices from the CDF file with their appropriate PMs (and
# hence assign the correct affinity).
 
  verbose && enter(verbose, "Identifying MM probes for these units");
# get mm indices
  isMm <- !isPm(this);
  verbose && exit(verbose);

# assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  indexPm <- sequenceInfo$y * dimension[1] + sequenceInfo$x + 1
  indexMmPutative <- indexPm + dimension[1]

# match up putative MM with actual list of MMs from CDF
  matches <- match(indexMmPutative, indices[isMm])
  indexMm <- indices[isMm][matches]
  notNA <- which(!is.na(indexMm))
  indexMm <- indexMm[notNA]

# now calculate affinities - code reused from compute.affinities() in
# gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  data(affinity.spline.coefs)

  affinity.basis.matrix <- ns(1:25, df = length(affinity.spline.coefs)/3)

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5])
  T13 <- 0
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10])
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15])

  apm <- vector("double", nbrOfSequences);
  amm <- vector("double", nbrOfSequences);

  if (verbose) {
    cat(verbose, "Progress: ");
    pb <- ProgressBar(stepLength=100/(nbrOfSequences/1000));
    reset(pb);
  }
  for (i in seq(along = apm)) {
    if (verbose && i %% 1000 == 0)
      increase(pb);
    # Get a 4x25 matrix with rows A, C, G, and T.
    charMtrx <- .Call("gcrma_getSeq", sequenceInfo$sequence[i], PACKAGE = "gcrma")

    A <- cbind(charMtrx[1, ] %*% affinity.basis.matrix, 
               charMtrx[2, ] %*% affinity.basis.matrix, 
               charMtrx[3, ] %*% affinity.basis.matrix)
    apm[i] <- A %*% affinity.spline.coefs
    if (charMtrx[1, 13] == 1) {
      amm[i] <- apm[i] + T13 - A13
    } else {
      if (charMtrx[4, 13] == 1) {
        amm[i] <- apm[i] + A13 - T13
      } else {
        if (charMtrx[3, 13]) {
          amm[i] <- apm[i] + C13 - G13
        } else {
          amm[i] <- apm[i] + G13 - C13
        }
      }
    }
  }
  verbose && cat(verbose, "");

  # create a vector to hold affinities and assign values to the appropriate
  # location in the vector
  affinities <- rep(NA, dimension[1]*dimension[2]);
  affinities[indexPm] <- apm;
  affinities[indexMm] <- amm[notNA];

  verbose && exit(verbose);

  # Saving to cache
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, affinities, comment=comment, dirs=dirs);

  verbose && exit(verbose);
  
  affinities;
}, private=TRUE)

############################################################################
# HISTORY:
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
