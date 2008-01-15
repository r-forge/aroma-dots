%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getBackgroundSD.GenePixData}
\alias{getBackgroundSD.GenePixData}
\alias{GenePixData.getBackgroundSD}
\alias{getBackgroundSD.GenePixData}
\alias{getBackgroundSD,GenePixData-method}


\title{Gets the standard deviation of the background pixels}

\usage{\method{getBackgroundSD}{GenePixData}(this, slides=NULL, channels=1:2, ...)}

\description{
 Gets the standard deviation of the background pixels.
}

\value{
  Returns a list of matrices that contain the estimated standard deviation
  of the pixels in the background region of the spots.
}

\details{
  The GenePix fields \code{"B532 SD"} (green) and \code{"B635 SD"} (red)
  are returned.
}

\examples{\dontrun{See help for getForegroundSD() for this class.}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
