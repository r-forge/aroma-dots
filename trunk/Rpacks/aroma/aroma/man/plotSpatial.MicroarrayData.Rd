%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.MicroarrayData}
\alias{plotSpatial.MicroarrayData}
\alias{MicroarrayData.plotSpatial}
\alias{plotSpatial.MicroarrayData}
\alias{plotSpatial,MicroarrayData-method}


\title{Plots a spatial representation of one field}

\description{
 Plots a spatial representation of one field.
 It is recommended to use the \code{plot} function instead of calling
 this method explicitly (see \code{\link[aroma:plot.MicroarrayData]{*plot}()}).
}

\usage{\method{plotSpatial}{MicroarrayData}(this, what, slide=1, include=NULL, exclude=NULL, col="auto", xlab=NULL, ylab="", axes=FALSE, xaxs="i", yaxs="i", pch="auto", grid=TRUE, log=NULL, ..., cex=NULL, style=NULL)}

\arguments{
 \item{what}{What to plot. Any field that can be retrieved by \code{extract},
  is accepted.}
 \item{slide}{The slide to be plotted.}
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded. If \code{\link[base]{NULL}} no spots are excluded.}
 \item{col}{The color(s) to be used for the plotted spots, i.e. for the
  spots \emph{after} inclusion and exclusion. If the value is
  \code{"redgreen"} a red to green palette is used.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
 \item{cex}{For internal use only! See above.}
 \item{pch}{For internal use only! See above.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  subplots(4)
  plotSpatial(ma)                   # Spatial plot of log ratios before.
  normalizeWithinSlide(ma, "p")     # Printtipwise lowess normalization.
  plotSpatial(ma)                   # Spatial plot of log ratios after.

  plotSpatial(ma, include=(abs(ma$M) > 2))
  points(ma, include=(abs(ma$M) > 2), col="red")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:plot.MicroarrayData]{*plot}()}.
  \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()}.
  \code{\link[aroma:highlight.MicroarrayData]{*highlight}()}.
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
