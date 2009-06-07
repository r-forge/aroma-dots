%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setSlideNames.MicroarrayData}
\alias{setSlideNames.MicroarrayData}
\alias{MicroarrayData.setSlideNames}
\alias{setSlideNames.MicroarrayData}
\alias{setSlideNames,MicroarrayData-method}


\title{Sets the names of the slides}

\usage{\method{setSlideNames}{MicroarrayData}(this, names, slides=NULL, ...)}

\description{
 Sets the names of the slides.
}

\arguments{
 \item{names}{A \code{\link[base]{vector}} of \code{\link[base]{character}} strings for the new slide names.}
 \item{slides}{A \code{\link[base]{vector}} (of the same length as \code{names}) of slides
    whose names should be set. If \code{\link[base]{NULL}}, all slides are considered.}
}

\value{
  Returns a \code{\link[base]{vector}} of the old slide names that were replaced.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getSlideNames.MicroarrayData]{*getSlideNames}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
