%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:02 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{keepSlides.MicroarrayData}
\alias{keepSlides.MicroarrayData}
\alias{MicroarrayData.keepSlides}
\alias{keepSlides.MicroarrayData}
\alias{keepSlides,MicroarrayData-method}


\title{Remove all but the specified slides}

\description{
  Remove all but the specified slides.

  Note that this method can also be used to reshuffle the slides in a
  certain order (and even duplicate some slides).
}

\usage{\method{keepSlides}{MicroarrayData}(this, slides, ...)}

\arguments{
  \item{slides}{\code{\link[base]{vector}} of slides to be kept. If \code{\link[base]{NULL}}, nothing is done.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:append.MicroarrayData]{*append}()} and \code{\link[aroma:removeSlides.MicroarrayData]{*removeSlides}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
