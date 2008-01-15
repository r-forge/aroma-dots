%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getArea.GenePixData}
\alias{getArea.GenePixData}
\alias{GenePixData.getArea}
\alias{getArea.GenePixData}
\alias{getArea,GenePixData-method}


\title{Gets the area of the spot}

\description{
 Gets the area of the spot. The area is the number of pixels identified (segmentated)
 to be within the spot limits.
}

\usage{\method{getArea}{GenePixData}(this, slides=NULL, include=NULL, ...)}

\arguments{
  \item{slides}{Specifying for which slides the spot area should
    be extracted. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{include}{The spots for which the spot area is returned.
      If \code{\link[base]{NULL}} all spots are considered.}
}

\value{
  Returns a \code{\link[base]{matrix}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
