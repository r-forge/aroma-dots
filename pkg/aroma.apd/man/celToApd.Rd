%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  celToApd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{celToApd}
\alias{celToApd.default}
\alias{celToApd}

\title{Generates an APD file from an Affymetrix CEL file}

\usage{
\method{celToApd}{default}(filename, apdFile=NULL, mapType="asChipType", writeMap=NULL, ..., verbose=FALSE)
}

\description{
  Generates an APD file from an Affymetrix CEL file.
}

\arguments{
  \item{filename}{The pathname of the CEL file.}
  \item{apdFile}{An optional pathname of the APD file, otherwise
    it will be the same as the CEL file with extension replaced
    with 'apd'.}
  \item{mapType}{The type of read map for the generated APD file.
    If \code{\link[base]{NULL}}, no remapping of the cell indices is done.
    If \code{"asChipType"}, the map type is the same as the chip type
    of the CEL file.
    If any other \code{\link[base]{character}} string, it sets the map type to that string.
    Note that there must be a APD map file with this type that can
    be found by \code{\link{findApdMap}}().
  }
  \item{writeMap}{An optional \emph{write} map \code{\link[base]{integer}} \code{\link[base]{vector}} used
    to remap the cell indices for optimal reading speed.  If \code{\link[base]{NULL}},
    the write map may be obtained from the read map file specified by
    the \code{mapType} argument.}
  \item{...}{Additional arguments passed to \code{\link{writeApd}}().}
  \item{verbose}{A \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns (invisibly) the pathname of the written APD file.
}

\examples{

library("R.utils")  ## Arguments

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Scan for existing CEL files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a) Scan for CEL files
files <- list.files(pattern="[.](cel|CEL)$")
files <- files[!file.info(files)$isdir]
if (length(files) > 0 && require("affxparser")) {
  cat("Create an optimal read map for CEL file:", files[1], "\n")
  cdffile <- findCdf(readCelHeader(files[1])$chiptype)
  res <- cdfToApdMap(cdffile)
 
  cat("Converting CEL file to APD file:", files[1], "\n")
  apdfile <- celToApd(files[1])
  cat("Created APD file:", apdfile, "\n")
  file.remove(apdfile)
 
  cat("Converting CEL file to APD file with an optimized read map:", files[1], "\n")
  apdfile <- celToApd(files[1], mapType=res$mapType)
  cat("Created APD file:", apdfile, "\n")
 
  writeMap <- invertMap(res$readMap)
  for (file in files[-1]) {
    cat("Converting CEL file to APD file with an optimized read map:", file, "\n")
    apdfile <- celToApd(file, mapType=res$mapType, writeMap=writeMap)
    cat("Created APD file:", apdfile, "\n")
  }
} # if (length(files) > 0)
}

\author{Henrik Bengtsson}

\seealso{
  To create an APD map file from a CDF file, see \code{\link{cdfToApdMap}}().
  To read an APD file, see \code{\link{readApd}}().
  To read an APD map file, see \code{\link{readApdMap}}().
}



\keyword{file}
\keyword{IO}
