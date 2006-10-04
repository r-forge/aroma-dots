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
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector of probe affinities, of length equal
#  to the total number of features on the array.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("computeAffinities", "AffymetrixCdfFile", function(this, paths=NULL, ..., verbose=FALSE) {

  chipType <- getChipType(this)

# try to find probe sequence file

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

  psFile <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
  if (is.null(psFile)) 
    throw("Could not locate probe sequence file: ", pattern);
  
# read in probe sequence data

  dimension <- getDimension(this)

  con <- file(psFile, "r")

# get header information and find columns containing (x,y) and sequence  
  hdr <- read.delim(file=con, header=FALSE, nrows=1, as.is=TRUE)
  hdr <- as.character(hdr)
  xcol <- grep("[Pp]robe [Xx]", hdr)
  ycol <- grep("[Pp]robe [Yy]", hdr)
  seqcol <- grep("[Pp]robe [Ss]equence", hdr)

  # check that file contains all three required columns: probe X and Y
  # and sequence
  if (is.null(xcol))
    throw("probe sequence file does not contain column named \"probe x\" (or any variation thereof): ", psFile);
  if (is.null(ycol))
    throw("probe sequence file does not contain column named \"probe y\" (or any variation thereof): ", psFile);
  if (is.null(seqcol))
    throw("probe sequence file does not contain column named \"probe sequence\" (or any variation thereof): ", psFile);
  

  # read everybody in using scan() with data types specified to speed things
  # up
  what <- vector("list", length(hdr))
  names(what) <- rep("dummy", length(hdr))
  what[xcol] <- list("numeric")
  what[ycol] <- list("numeric")
  what[seqcol] <- list("character")
  names(what)[xcol] <- "x"
  names(what)[ycol] <- "y"
  names(what)[seqcol] <- "sequence"
  sep <- "\t"

  sequenceInfo <- scan(file=con, what=what, sep=sep, quiet=TRUE) 
  sequenceInfo$x <- as.numeric(sequenceInfo$x)
  sequenceInfo$y <- as.numeric(sequenceInfo$y)
  indexPm <- sequenceInfo$y * dimension[1] + sequenceInfo$x + 1
  
  close(con)

  indices <- unlist(getCellIndices(this))

# probe sequence tab-separated files generally only contain the PM
# sequence (since we know that MM sequence always has base 13 complemented).
# This is true for expression arrays, but possibly not for genotyping arrays.
# We will calculate affinity for PM and MM anyway, and then
# try to match MM indices from the CDF file with their appropriate PMs (and
# hence assign the correct affinity).
  
# get mm indices
  isMm <- !isPm(this)

# assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  indexMmPutative <- indexPm+dimension[1]

# match up putative MM with actual list of MMs from CDF
  matches <- match(indexMmPutative, indices[isMm])
  indexMm <- indices[isMm][matches]
  notNA <- which(!is.na(indexMm))
  indexMm <- indexMm[notNA]

# now calculate affinities - code reused from compute.affinities() in
# gcrma

  require(gcrma, quietly = TRUE)
  require(splines, quietly = TRUE)
  require(matchprobes, quietly = TRUE)
  data(affinity.spline.coefs)

  affinity.basis.matrix <- ns(1:25, df = length(affinity.spline.coefs)/3)

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5])
  T13 <- 0
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10])
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15])

  apm <- vector("numeric", length(sequenceInfo$sequence));
  amm <- vector("numeric", length(sequenceInfo$sequence));

  for (i in seq(along = apm)) {
    charMtrx <- .Call("gcrma_getSeq", sequenceInfo$sequence[i], PACKAGE = "gcrma")
    A <- cbind(charMtrx[1, ] %*% affinity.basis.matrix, charMtrx[2,
                                                                 ] %*% affinity.basis.matrix, charMtrx[3, ] %*% affinity.basis.matrix)
    apm[i] <- A %*% affinity.spline.coefs
    if (charMtrx[1, 13] == 1) {
      amm[i] <- apm[i] + T13 - A13
    }
    else {
      if (charMtrx[4, 13] == 1) {
        amm[i] <- apm[i] + A13 - T13
      }
      else {
        if (charMtrx[3, 13]) {
          amm[i] <- apm[i] + C13 - G13
        }
        else {
          amm[i] <- apm[i] + G13 - C13
        }
      }
    }
  }


  # create a vector to hold affinities and assign values to the appropriate
  # location in the vector
  affinities <- vector("double", dimension[1]*dimension[2]);
  
  affinities[indexPm] <- apm;
  affinities[indexMm] <- amm[notNA];
  
  return(affinities);
  
})

############################################################################
# HISTORY:
# 2006-10-04
# o Debugged, tested, docs/comments added.
# 2006-10-01
# o Created.
############################################################################
