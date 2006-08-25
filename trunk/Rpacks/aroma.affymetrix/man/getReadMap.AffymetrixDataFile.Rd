%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.MAP.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getReadMap.AffymetrixDataFile}
\alias{getReadMap.AffymetrixDataFile}
\alias{AffymetrixDataFile.getReadMap}
\alias{getReadMap.AffymetrixDataFile}
\alias{getReadMap,AffymetrixDataFile-method}


\title{Gets the probe read map for a data file}

\description{
  Gets the probe read map for a data file.
}

\usage{\method{getReadMap}{AffymetrixDataFile}(this, inverse=FALSE, ...)}

\arguments{
 \item{inverse}{If \code{\link[base:logical]{FALSE}}, the read map is returned, otherwise the write map.}
 \item{...}{Not used.}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:setApdMap.AffymetrixDataFile]{*setApdMap}()}.
  \code{\link[aroma.affymetrix:hasApdMap.AffymetrixDataFile]{*hasApdMap}()}.
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
