%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.UnitModel}
\alias{readUnits.UnitModel}
\alias{UnitModel.readUnits}
\alias{readUnits.UnitModel}
\alias{readUnits,UnitModel-method}


\title{Reads data unit by unit}

\description{
 Reads data unit by unit for all or a subset of units (probesets).
}

\usage{\method{readUnits}{UnitModel}(this, units=NULL, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Arguments passed to \code{getCellIndices()} of the
    \code{\link{AffymetrixCdfFile}} class (if \code{cdf} was not specified),
    as well as \code{getUnitIntensities()} of the
    \code{\link{AffymetrixCelSet}} of the input data set.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached cell indices as well as data is re-read.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{readUnits()} of
 \code{\link{AffymetrixCelSet}} returns.
}

\details{
  The output structure is shaped by the \code{\link[base]{list}} structure returned
  by \code{\link[aroma.affymetrix:getCellIndices.UnitModel]{*getCellIndices}()}.  The default is to return whatever
  the CDF returns, but by overriding the latter, what cells are read
  and what group-structure each unit has can be fully customized
  by any subclass.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getCellIndices.UnitModel]{*getCellIndices}()}.
  For more information see \code{\link{UnitModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
