%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  updateApdUnits.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{updateApdUnits}
\alias{updateApdUnits.default}
\alias{updateApdUnits}

\title{Updates an Affymetrix probe data (APD) file by units (probesets)}

\usage{
\method{updateApdUnits}{default}(filename, units=NULL, data, ..., cdf=NULL, stratifyBy=c("nothing", "pmmm",
  "pm", "mm"), verbose=FALSE)
}

\description{
  Updates an Affymetrix probe data (APD) file by units (probesets) by using the unit and group definitions in the
  corresponding Affymetrix CDF file.
}

\arguments{
  \item{filename}{The filename of the APD file.}
  \item{units}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of unit indices specifying which
    units to be read.  If \code{\link[base]{NULL}}, all units are updated.}
  \item{data}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of data elements to be assigned.}
  \item{...}{Additional arguments passed to \code{\link{updateApd}}(),
    e.g. \code{writeMap}.}
  \item{cdf}{A \code{\link[base]{character}} filename of a CDF file, or a CDF \code{\link[base]{list}}
    structure.  If \code{\link[base]{NULL}}, the CDF file is searched for by
    \code{\link[affxparser]{findCdf}}.}
  \item{stratifyBy}{Argument passed to low-level method
    \code{\link[affxparser]{readCdfCellIndices}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson}

\examples{\dontrun{#See ?createApd for an example.}}

\seealso{
  \code{\link{readApdUnits}}() to read unit by units.
  \code{\link{updateApd}}() to update cell by cell.
}



\keyword{file}
\keyword{IO}
