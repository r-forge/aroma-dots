%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotfinderData.R
% 
% on Tue Jan 15 18:37:29 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.SpotfinderData}
\alias{plotSpatial.SpotfinderData}
\alias{SpotfinderData.plotSpatial}
\alias{plotSpatial.SpotfinderData}
\alias{plotSpatial,SpotfinderData-method}


\title{Creates a spatial plot of a slide}

\description{
 Creates a spatial plot of a slide.
 Note that older versions of the
 Spot software did not generate result file containing information about
 the spatial coordinates of the spots. If this is the case, a standard
 spatial plot is generated.
}

\usage{\method{plotSpatial}{SpotfinderData}(this, what="meanRatio", slide=1, pch=20, yaxt=NULL, ...)}

\arguments{
  \item{slide}{The slide to be plotted.}
  \item{pch}{The spot symbol. Default is \code{20} (solid disk).}
  \item{yaxt, ...}{Parameters accepted by \code{plot}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{# For an example see help(SpotfinderData).}}

\seealso{
  For more information see \code{\link{SpotfinderData}}.
}
\keyword{internal}
\keyword{methods}
