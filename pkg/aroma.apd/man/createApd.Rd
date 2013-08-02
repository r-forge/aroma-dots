%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  createApd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{createApd}
\alias{createApd.default}
\alias{createApd}

\title{Creates an Affymetrix probe data (APD) file}

\usage{
\method{createApd}{default}(filename, nbrOfCells, dataType=c("float", "double", "integer", "short",
  "byte"), chipType=NULL, mapType=NULL, ..., verbose=FALSE, .checkArgs=TRUE)
}

\description{
  Creates an Affymetrix probe data (APD) file.

  An Affymetrix probe data (APD) structure can hold a header and
  a numeric data vector.  Since the APD structure is kept on file
  all the time, the number of elements in the data vector is only
  limited by the file system and not the amount of system memory
  available.  For more details, see the \code{\link[R.huge]{FileVector}}
  class (and its superclass), which is used internally.
}

\arguments{
  \item{filename}{The filename of the APD file.}
  \item{nbrOfCells}{The number of cells (probes) data elements the
     APD file structure should hold.}
  \item{dataType}{The data type of the data elements.}
  \item{chipType}{An (optional) \code{\link[base]{character}} string specifying the
     chip type.}
  \item{mapType}{An (optional) \code{\link[base]{character}} string specifying the
     probe-index map.  Use by \code{\link{findApdMap}}() to find read map.}
  \item{...}{Additional named arguments added to the header of the
     APD file structure.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.checkArgs}{If \code{\link[base:logical]{TRUE}}, arguments are checked, otherwise not.}
}

\value{
  Returns (invisibly) the pathname of the file created.
}

\section{Data type}{
  Valid data types are: byte (1 byte), short (2 bytes), integer (4 bytes),
  float (4 bytes), and double (8 bytes).

  Note that in Affymetrix CEL files, the probe intensities as well as
  the standard deviations are stored as floats (4 bytes) and not doubles
  (8 bytes).  This is why, the default data type is \code{"float"}.
}

\author{Henrik Bengtsson}

\examples{

for (zzz in 0) {

# Float precision
.Machine$float.eps <- (2^((8-4)*8)*.Machine$double.eps)
tol <- .Machine$float.eps ^ 0.5

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Create an Affymetrix Probe Signal (APD) file for a
#    'Mapping50K_Hind240' with 1600-by-1600 probes.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipType <- "Mapping50K_Hind240"
nbrOfCells <- 1600^2

pathname <- paste(tempfile(), "apd", sep=".")
createApd(pathname, nbrOfCells=nbrOfCells, chipType=chipType)

# File size
cat("File name:", pathname, "\n")
cat("File size:", file.info(pathname)$size, "bytes\n")
cat("APD header:\n")
header <- readApdHeader(pathname)
print(header)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Update the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cells <- c(1, 1100:1120, 220:201, 998300:999302)
signals <- log(cells+1, base=2)  # Fake signals
updateApd(pathname, indices=cells, data=signals)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Re-read the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
apd <- readApd(pathname, indices=cells)

# Signals in APD files are stored as floats (since this is
# the precision in CEL files).
stopifnot(all.equal(signals, apd$intensities, tolerance=tol))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Re-read the signals for a subset of probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdfFile <- findCdf(chipType)
if (length(cdfFile) == 0)
  break

apd <- readApdUnits(pathname, units=100:104)

# Sample new data (with one decimal precision)
apd2 <- lapply(apd, function(unit) {
  lapply(unit, function(groups) {
    n <- length(groups$intensities)
    values <- as.integer(runif(n, max=655350))/10
    list(intensities=values)
  })
})

# Update APD file with new data
updateApdUnits(pathname, units=100:104, data=apd2)

# Re-read data to verify correctness
apd <- readApdUnits(pathname, units=100:104)

# Signals in APD files are stored as floats (since this is
# the precision in CEL files).
stopifnot(all.equal(apd2, apd, tolerance=tol))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Clean up
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
file.remove(pathname)

rm(apd2)

} # for (zzz in 0)

rm(zzz, tol, chipType, nbrOfCells, pathname, header, cells, signals,
                                                        apd, cdfFile)
}

\seealso{
  \code{\link{updateApd}}() and \code{\link{readApd}}().
  To find a map of a certain type, see \code{\link{findApdMap}}().
}



\keyword{file}
\keyword{IO}
