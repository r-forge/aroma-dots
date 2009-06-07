%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getIndex.Layout}
\alias{getIndex.Layout}
\alias{Layout.getIndex}
\alias{getIndex.Layout}
\alias{getIndex,Layout-method}


\title{Gets the index of a spot given its location}

\usage{\method{getIndex}{Layout}(this, gridRow, gridColumn=NULL, spotRow=NULL, spotColumn=NULL, ...)}

\description{
 Gets the index of a spot given its location. The location can either be
 a vector containing the grid row and the grid column and the spot row
 and the spot column in that grid, or it can be the same fields as
 seperate arguments.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # Example 1
  layout <- Layout(4,4, 18,18)
  idx <- getIndex(layout, 2, 3, 4, 3)     # 2001
  idx <- getIndex(layout, c(2, 3, 4, 3))  # 2001 (equivalent)
  loc <- getLocation(layout, idx)         # 2 3 4 3

  # Example 2
  SMA$loadData(c("mouse.data", "mouse.setup"))
  raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
  ma <- getSignal(raw)
  layout <- getLayout(ma)

  plotSpatial(ma)

  # Highlights spot number 2462
  idx <- 2462
  highlight(ma, idx, col="purple")

  # Highlights the spot at grid (2,3) and its spot (4,3)
  idx <- getIndex(layout, 2, 3, 4, 3);  # Spot #2460
  highlight(ma, idx, col="purple")

  # Highlights all spots in grid (1,2)
  idx <- getIndices(layout, 1,2, NULL,NULL)
  highlight(ma, idx, col="purple")
}

\seealso{
  This method corresponds to image2id (see \code{\link[sma]{id2image}}) in
  the sma package.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
