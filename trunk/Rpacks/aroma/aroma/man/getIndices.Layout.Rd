%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getIndices.Layout}
\alias{getIndices.Layout}
\alias{Layout.getIndices}
\alias{getIndices.Layout}
\alias{getIndices,Layout-method}


\title{Gets the indices of the spots at the given locations}

\usage{\method{getIndices}{Layout}(this, gridRows=NULL, gridColumns=NULL, spotRows=NULL, spotColumns=NULL, ...)}

\arguments{
 \item{gridRows}{The grid rows to be included. If \code{\link[base]{NULL}} all are included.}
 \item{gridColumns}{The grid column to be included. If \code{\link[base]{NULL}} all are included.}
 \item{spotRows}{The spot rows to be included. If \code{\link[base]{NULL}} all are included.}
 \item{spotColumns}{The spot column to be included. If \code{\link[base]{NULL}} all are included.}
}

\description{
 Gets the indices of the spots at the given locations.
 This method is very useful for instance when one would like to look for
 spatial effects such as edge effects etc.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.data", "mouse.setup"))
  raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
  ma <- getSignal(raw)
  layout <- getLayout(ma)

  plotSpatial(ma)

  # Highlights all spots in grid (1,2)
  idx <- getIndices(layout, gridRows=1, gridColumns=2);
  highlight(ma, idx, col="purple")

  # Highlight all spots in column 1, 6 and 12 in
  # grid (3,2), (3,3), (4,2) and (4,3).
  idx <- getIndices(layout, gridRows=3:4, gridColumns=2:3, spotColumns=c(1,6,12));
  highlight(ma, idx, col="orange")

  # Highlight all "alley spots" of each printtip group, i.e. those spots that
  # do *not* have eight neighbors and do a background "alley" next to them.
  alley <- getIndices(layout, spotRows=c(1,layout$nspot.r))
  alley <- union(alley, getIndices(layout, spotColumns=c(1,layout$nspot.c)));
  highlight(ma, alley, col="pink")
}
\keyword{internal}
\keyword{methods}
