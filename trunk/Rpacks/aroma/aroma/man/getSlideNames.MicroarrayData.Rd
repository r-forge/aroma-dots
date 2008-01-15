%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSlideNames.MicroarrayData}
\alias{getSlideNames.MicroarrayData}
\alias{MicroarrayData.getSlideNames}
\alias{getSlideNames.MicroarrayData}
\alias{getSlideNames,MicroarrayData-method}


\title{Gets the names of the slides}

\usage{\method{getSlideNames}{MicroarrayData}(this, slides=NULL, ...)}

\description{
 Gets the names of the slides.
}

\arguments{
 \item{slides}{A \code{\link[base]{vector}} of slides whose names should be returned.
    If \code{\link[base]{NULL}}, all slides are considered.}
}

\value{
  Returns a \code{\link[base]{vector}} of character strings or \code{\link[base]{NULL}}, if no slide names are
  set.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:setSlideNames.MicroarrayData]{*setSlideNames}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
