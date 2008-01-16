%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% on Tue Jan 15 18:36:51 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSpotPosition.ScanAlyzeData}
\alias{getSpotPosition.ScanAlyzeData}
\alias{ScanAlyzeData.getSpotPosition}
\alias{getSpotPosition.ScanAlyzeData}
\alias{getSpotPosition,ScanAlyzeData-method}


\title{Gets physical positions of the spots}

\description{
 Gets physical positions (in pixels) of the spots on one or several
 slides.
}

\usage{\method{getSpotPosition}{ScanAlyzeData}(this, slides=NULL, index=NULL, ...)}

\arguments{
  \item{slides}{Specifying which for slides the spot positions should
    be extracted. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{index}{The spots for which the position is returned.
    If \code{\link[base]{NULL}} all spots are considered.}
}

\value{Returns a \code{\link{SpotPosition}} object containing the
  positions of the spots on the specified slides.}

\details{
  The ScanAlyze software does not return the center position of a spot,
  but the \code{TOP} \code{LEFT} \code{BOT} and \code{RIGHT} coordinates.
  This method assumes that the center of the spot is in the center of
  this box.
}

\examples{
  sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))

  # Gets the spot positions
  xy <- getSpotPosition(sa)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ScanAlyzeData}}.
}
\keyword{internal}
\keyword{methods}
