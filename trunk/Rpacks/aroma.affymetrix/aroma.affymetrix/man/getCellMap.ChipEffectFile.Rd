%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellMap.ChipEffectFile}
\alias{getCellMap.ChipEffectFile}
\alias{ChipEffectFile.getCellMap}
\alias{getCellMap.ChipEffectFile}
\alias{getCellMap,ChipEffectFile-method}


\title{Gets a (unit, group, cell) index map}

\description{
 Gets a (unit, group, cell) index map.
}

\usage{\method{getCellMap}{ChipEffectFile}(this, units=NULL, force=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units for which the map should be returned.
     If \code{\link[base]{NULL}}, all units are considered.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached cell indices are ignored.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{data.frame}} with \code{\link[base]{integer}} columns \code{unit}, \code{group},
 and \code{cell}.
}

\examples{\dontrun{
     unit group cell
  1    20     1   20
  2    21     1   21
  3    22     1   22
  4    23     1   23
  5    24     1   24
  6    25     1   25
  7    25     2   27
  8    26     1   29
  9    26     2   31
}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
