%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSpotPosition.GenePixData}
\alias{getSpotPosition.GenePixData}
\alias{GenePixData.getSpotPosition}
\alias{getSpotPosition.GenePixData}
\alias{getSpotPosition,GenePixData-method}


\title{Gets physical positions of the spots}

\description{
 Gets physical positions (in pixels) of the spots on one or several
 slides.
}

\usage{\method{getSpotPosition}{GenePixData}(this, slides=NULL, index=NULL, ...)}

\arguments{
  \item{slides}{Specifying for which slides the spot positions should
    be extracted. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{index}{The spots for which the position is returned.
    If \code{\link[base]{NULL}} all spots are considered.}
}

\value{Returns a \code{\link{SpotPosition}} object containing the
  positions of the spots on the specified slides.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Gets the spot positions
  xy <- getSpotPosition(gpr)
}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
