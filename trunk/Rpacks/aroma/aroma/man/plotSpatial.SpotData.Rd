%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% on Tue Jan 15 18:37:20 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.SpotData}
\alias{plotSpatial.SpotData}
\alias{SpotData.plotSpatial}
\alias{plotSpatial.SpotData}
\alias{plotSpatial,SpotData-method}


\title{Creates a spatial plot of a slide}

\description{
  Creates a spatial plot of a slide.
 Note that older versions of the
 Spot software did not generate result file containing information about
 the spatial coordinates of the spots. If this is the case, a standard
 spatial plot is generated.
}

\usage{\method{plotSpatial}{SpotData}(this, what=NULL, slide=1, include=NULL, exclude=NULL, pch=20, yaxt=NULL, cex=NULL, col="auto", palette="redgreen", A.range=c(0, 16), M.range=c(-1, 1), look=c("real", "classic"), style=NULL, ...)}

\arguments{
  \item{slide}{The slide to be plotted.}
  \item{pch}{The spot symbol. Default is \code{20} (solid disk).}
  \item{yaxt, ...}{Parameters accepted by \code{plot}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{For an example see help(SpotData).}}

\seealso{
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
