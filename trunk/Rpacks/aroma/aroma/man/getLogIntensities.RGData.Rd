%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLogIntensities.RGData}
\alias{getLogIntensities.RGData}
\alias{RGData.getLogIntensities}
\alias{getLogIntensities.RGData}
\alias{getLogIntensities,RGData-method}

\alias{getA}

\title{Calculates the log-intensitites (A values)}

\description{
  Calculates the log-intensitites (A values) for each specified array.
  Log-base 2 (two) is used.
}

\usage{\method{getLogIntensities}{RGData}(this, slides=NULL, ...)}

\arguments{
  \item{slides}{A \code{\link[base]{vector}} of \code{\link[base]{integer}}s indicating which slides to be
    considered. If \code{\link[base]{NULL}}, all slides are considered.}
}

\value{
  Returns a NxK \code{\link[base]{matrix}} where N is the number of spots and K is the
  number of (specified) arrays.
}

\seealso{
  \code{\link[aroma:getLogRatios.RGData]{*getLogRatios}()}.
  For more information see \code{\link{RGData}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
