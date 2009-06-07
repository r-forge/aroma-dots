%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{keepSpots.MicroarrayData}
\alias{keepSpots.MicroarrayData}
\alias{MicroarrayData.keepSpots}
\alias{keepSpots.MicroarrayData}
\alias{keepSpots,MicroarrayData-method}


\title{Remove all but the specified spots from all slides}

\description{
  Remove all but the specified spots from all slides.

  Note that this method can also be used to reshuffle the spots in a
  certain order (and even duplicate some spots).
}

\usage{\method{keepSpots}{MicroarrayData}(this, index, ...)}

\arguments{
  \item{index}{\code{\link[base]{vector}} of spot indices to be kept. If \code{\link[base]{NULL}}, nothing
   is done.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:removeSpots.MicroarrayData]{*removeSpots}()}.
  For similar operations on slides, see \code{\link[aroma:keepSlides.MicroarrayData]{*keepSlides}()} and
  \code{\link[aroma:removeSlides.MicroarrayData]{*removeSlides}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
