%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{removeSpots.MicroarrayData}
\alias{removeSpots.MicroarrayData}
\alias{MicroarrayData.removeSpots}
\alias{removeSpots.MicroarrayData}
\alias{removeSpots,MicroarrayData-method}


\title{Remove specified spots}

\description{
  Remove specified spots.
}

\usage{\method{removeSpots}{MicroarrayData}(this, index, ...)}

\arguments{
 \item{spots}{\code{\link[base]{vector}} of spots to be removed.
              If \code{\link[base]{NULL}}, all spots are removed.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:keepSpots.MicroarrayData]{*keepSpots}()}.
  For similar operations on slides, see \code{\link[aroma:keepSlides.MicroarrayData]{*keepSlides}()} and
  \code{\link[aroma:removeSlides.MicroarrayData]{*removeSlides}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
