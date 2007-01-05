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
 Reads data unit by unit for all or a subset of units (probeset)
 specially structured for this PLM.
}

\usage{\method{readUnits}{UnitModel}(this, units=NULL, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units to be read. If \code{\link[base]{NULL}}, all units are read.}
  \item{...}{Arguments passed to \code{getCellIndices()} of the
    \code{\link{AffymetrixCdfFile}} class (if \code{cdf} was not specified),
    but also to the \code{readUnits()} method of the
    \code{\link{AffymetrixCelSet}} class.}
}

\value{
 Returns the \code{\link[base]{list}} structure that \code{readUnits()} of
 \code{\link{AffymetrixCelSet}} returns.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:updateUnits.UnitModel]{*updateUnits}()}.
  For more information see \code{\link{UnitModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
