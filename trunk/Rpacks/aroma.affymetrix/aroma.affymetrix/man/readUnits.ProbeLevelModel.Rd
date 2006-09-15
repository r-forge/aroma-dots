%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.ProbeLevelModel}
\alias{readUnits.ProbeLevelModel}
\alias{ProbeLevelModel.readUnits}
\alias{readUnits.ProbeLevelModel}
\alias{readUnits,ProbeLevelModel-method}


\title{Reads data unit by unit}

\description{
 Reads data unit by unit for all or a subset of units (probeset)
 specially structured for this PLM.
}

\usage{\method{readUnits}{ProbeLevelModel}(this, units=NULL, ..., verbose=FALSE)}

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
  \code{\link[aroma.affymetrix:updateUnits.ProbeLevelModel]{*updateUnits}()}.
  For more information see \code{\link{ProbeLevelModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
