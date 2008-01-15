%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLogRatios.RGData}
\alias{getLogRatios.RGData}
\alias{RGData.getLogRatios}
\alias{getLogRatios.RGData}
\alias{getLogRatios,RGData-method}

\alias{getM}

\title{Calculates the log-ratios (M values)}

\description{
  Calculates the log-ratios (M values) for each specified array.
  Log-base 2 (two) is used.
}

\usage{\method{getLogRatios}{RGData}(this, slides=NULL, ...)}

\arguments{
  \item{slides}{A \code{\link[base]{vector}} of \code{\link[base]{integer}}s indicating which slides to be
    considered. If \code{\link[base]{NULL}}, all slides are considered.}
}

\value{
  Returns a NxK \code{\link[base]{matrix}} where N is the number of spots and K is the
  number of (specified) arrays.
}

\seealso{
  \code{\link[aroma:getLogIntensities.RGData]{*getLogIntensities}()}.
  For more information see \code{\link{RGData}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
