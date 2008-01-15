%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{highlight.MicroarrayData}
\alias{highlight.MicroarrayData}
\alias{MicroarrayData.highlight}
\alias{highlight.MicroarrayData}
\alias{highlight,MicroarrayData-method}


\title{Highlight spots in last plot}

\description{
 Highlight spots in last plot.
}

\usage{\method{highlight}{MicroarrayData}(this, include=NULL, exclude=NULL, col="black", ...)}

\arguments{
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded.
  If \code{\link[base]{NULL}} no spots are excluded.}
 \item{pch}{The point style to be used. Default value is 176 (a small open circle). For more information see \code{\link[graphics]{par}}.}
 \item{cex}{The scale factor to be used. If \code{\link[base]{NULL}} the system default scale factor will be used. For more information see \code{\link[graphics]{par}}.}
 \item{col}{The color(s) to be used for the plotted spots, i.e. for the
  spots \emph{after} inclusion and exclusion.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  # Get all spots in the grid at row 1 and at column 2.
  layout <- getLayout(ma)
  idx <- getIndices(layout, 1,2, NULL,NULL)

  # Plot the data both highlighted and non-highlighted.
  layout(matrix(1:4, ncol=2, byrow=TRUE))
  plot(ma)
  plotSpatial(ma)
  plot(ma)
  highlight(ma, idx, col="purple")  # Highlight all spots in grid (1,2).
  plotSpatial(ma)
  highlight(ma, idx, col="purple")  # Highlight all spots in grid (1,2).
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  \code{\link[graphics]{points}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
