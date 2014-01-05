%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readApdUnits.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readApdUnits}
\alias{readApdUnits.default}
\alias{readApdUnits}

\title{Reads Affymetrix probe data (APD) as units (probesets)}

\usage{
\method{readApdUnits}{default}(filenames, units=NULL, ..., transforms=NULL, cdf=NULL,
  stratifyBy=c("nothing", "pmmm", "pm", "mm"), addDimnames=FALSE, readMap="byMapType",
  dropArrayDim=TRUE, verbose=FALSE)
}

\description{
  Reads Affymetrix probe data (APD) as units (probesets) by using the unit and group definitions in the
  corresponding Affymetrix CDF file.

  If more than one APD file is read, all files are assumed to be of
  the same chip type, and have the same read map, if any.
  It is not possible to read APD files of different types at the
  same time.
}

\arguments{
  \item{filenames}{The filenames of the APD files.  All APD files must
    be of the same chip type.}
  \item{units}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of unit indices specifying which units
    to be read.  If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Additional arguments passed to \code{\link{readApd}}().}
  \item{transforms}{A \code{\link[base]{list}} of exactly \code{length(filenames)}
    \code{\link[base]{function}}s.  If \code{\link[base]{NULL}}, no transformation is performed.
    Values read are passed through the corresponding transform
    function before being returned.}
  \item{cdf}{A \code{\link[base]{character}} filename of a CDF file, or a CDF \code{\link[base]{list}}
    structure.  If \code{\link[base]{NULL}}, the CDF file is searched for by
    \code{\link[affxparser]{findCdf}} first starting from the current directory
    and then from the directory where the first APD file is.}
  \item{stratifyBy}{Argument passed to low-level method
    \code{\link[affxparser]{readCdfCellIndices}}.}
  \item{addDimnames}{If \code{\link[base:logical]{TRUE}}, dimension names are added to arrays,
    otherwise not.  The size of the returned APD structure in bytes
    increases by 30-40\% with dimension names.}
  \item{readMap}{A \code{\link[base]{vector}} remapping cell indices to file indices.
    If \code{"byMapType"}, the read map of type according to APD header
    will be search for and read.  It is much faster to specify the
    read map explicitly compared with searching for it each time.
    If \code{\link[base]{NULL}}, no map is used.}
  \item{dropArrayDim}{If \code{\link[base:logical]{TRUE}} and only one array is read, the elements of
    the group field do \emph{not} have an array dimension.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  A named \code{\link[base]{list}} where the names corresponds to the names of the units
  read.  Each element of the \code{\link[base]{list}} is in turn a \code{\link[base]{list}} structure
  with groups (aka blocks).
}

\section{Speed}{
  Since the cell indices are semi-randomized across the array and
  with units (probesets), it is very unlikely that the read will
  consist of subsequent cells (which would be faster to read).
  However, the speed of this method, which uses \code{\link[R.huge]{FileVector}}
  to read data, is comparable to the speed of
  \code{\link[affxparser]{readCelUnits}}, which uses the Fusion SDK
  (\code{\link[affxparser]{readCel}}) to read data.
}

\author{Henrik Bengtsson}

\examples{

library("R.utils") # Arguments

verbose <- Arguments$getVerbose(TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Scan for existing CEL files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a) Scan current directory for CEL files
files <- list.files(pattern="[.](cel|CEL)$")
files <- files[!file.info(files)$isdir]

if (length(files) > 0 && require("affxparser")) {
  # b) Corresponding APD filenames
  celNames <- files
  apdNames <- gsub(pattern, ".apd", files)
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Copy the probe intensities from a CEL to an APD file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in 1) {
    verbose && enter(verbose, "Reading CEL file #", kk)
    cel <- readCel(celNames[kk])
    verbose && exit(verbose)
 
    if (!file.exists(apdNames[kk])) {
      verbose && enter(verbose, "Creating APD file #", kk)
      chipType <- cel$header$chiptype
      writeApd(apdNames[kk], data=cel$intensities, chipType=chipType)
      verbose && exit(verbose)
    }
 
    verbose && enter(verbose, "Verifying APD file #", kk)
    apd <- readApd(apdNames[kk])
    verbose && exit(verbose)
    stopifnot(identical(apd$intensities, cel$intensities))
  }
 
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Read a subset of the units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  units <- c(1, 20:205)
  cel <- readCelUnits(celNames[1], units=units)
  apd <- readApdUnits(apdNames[1], units=units)
  stopifnot(identical(apd, cel))
 
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. The same, but stratified on PMs and MMs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  apd <- readApdUnits(apdNames[1], units=units, stratifyBy="pmmm",
                                                addDimnames=TRUE)
} # if (length(files) > 0)
}

\seealso{
  To read CEL units, \code{\link[affxparser]{readCelUnits}}.
  Internally, the \code{\link{readApd}}() method is used for read probe data,
  and \code{\link{readApdMap}}(), if APD file has a map type specified and
  the read map was not given explicitly.
}

\references{
  [1] Affymetrix Inc, Affymetrix GCOS 1.x compatible file formats,
      June 14, 2005.
      \url{http://www.affymetrix.com/support/developer/}
}



\keyword{file}
\keyword{IO}
