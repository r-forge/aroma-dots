%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% on Tue Jan 15 18:37:20 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSpotPosition.SpotData}
\alias{getSpotPosition.SpotData}
\alias{SpotData.getSpotPosition}
\alias{getSpotPosition.SpotData}
\alias{getSpotPosition,SpotData-method}


\title{Gets physical positions of the spots}

\description{
 Gets physical positions (in pixels) of the spots on one or several
 slides.
}

\usage{\method{getSpotPosition}{SpotData}(this, slides=NULL, index=NULL, ...)}

\arguments{
  \item{slides}{Specifying which for slides the spot positions should
    be extracted. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{index}{The spots for which the position is returned.
    If \code{\link[base]{NULL}} all spots are considered.}
}

\value{Returns a \code{\link{SpotPosition}} object containing the
  positions of the spots on the specified slides.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
