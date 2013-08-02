%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readApd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readApd}
\alias{readApd.default}
\alias{readApd}

\title{Reads an Affymetrix probe data (APD) file}

\usage{
\method{readApd}{default}(filename, indices=NULL, readMap="byMapType", name=NULL, ..., verbose=FALSE,
  .checkArgs=TRUE)
}

\description{
  Reads an Affymetrix probe data (APD) file.
}

\arguments{
  \item{filename}{The filename of the APD file.}
  \item{indices}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of cell (probe) indices
    specifying what cells to read.  If \code{\link[base]{NULL}}, all are read.}
  \item{readMap}{A \code{\link[base]{vector}} remapping cell indices to file indices.
    If \code{"byMapType"}, the read map of type according to APD header
    will be search for and read.  It is much faster to specify the
    read map explicitly compared with searching for it each time.
    If \code{\link[base]{NULL}}, no map is used.}
  \item{name}{The name of the data field.
    If \code{\link[base]{NULL}}, the APD header \code{name} is used.  If not specified,
    it defaults to \code{"intensities"}.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.checkArgs}{If \code{\link[base:logical]{TRUE}}, arguments are validated, otherwise not.}
}

\value{
  A named \code{\link[base]{list}} with the two elements \code{header} and
  \code{data}.  The header is in turn a \code{\link[base]{list}} structure and
  the second is a \code{\link[base]{numeric}} \code{\link[base]{vector}} holding the queried data.
}

\details{
  To read one \emph{large} contiguous block of elements is faster than
  to read individual elements one by one.  For this reason, internally
  more elements than requested may be read and therefore allocation more
  memory than necessary.  This means, in worst case \eqn{N} elements
  may read allocation \eqn{N*8} bytes of \R memory, although only two
  elements are queried.  However, to date even with the largest arrays
  from Affymetrix this will still only require tens of megabytes of
  \emph{temporary} memory.  For instance, Affymetrix Mapping 100K arrays
  holds 2,560,000 probes requiring 20Mb of temporary memory.
}

\section{Remapping indices}{
  Argument \code{readMap} can be used to remap indices.  For instance,
  the indices of the probes can be reorder such that the probes within
  a probeset is in a contiguous set of probe indices.  Then, given that
  the values are stored in such an order, when reading complete probesets,
  data will be access much faster from file than if the values were
  scatter all over the file.

  Example of speed improvements.  Reading all 40000 values in units
  1001 to 2000 of an Affymetrix Mapping 100K Xba chip is more than
  10-30 times faster with mapping compared to without.
}

\section{File format}{
  The file format of an APD file is identical to the file format of an
  \code{\link[R.huge]{FileVector}}.
}

\author{Henrik Bengtsson}

\examples{\dontrun{#See ?createApd for an example.}}

\seealso{
  \code{\link{createApd}}() and \code{\link{updateApd}}().
  See also \code{\link{readApdHeader}}().
  To create a cell-index read map from an CDF file, see
  \code{\link[affxparser]{readCdfUnitsWriteMap}}.
}



\keyword{file}
\keyword{IO}
