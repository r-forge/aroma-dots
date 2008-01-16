%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plot.MicroarrayData}
\alias{plot.MicroarrayData}
\alias{MicroarrayData.plot}
\alias{plot.MicroarrayData}
\alias{plot,MicroarrayData-method}


\title{Plots spatial representation of a microarray spot statistics}

\description{
 Plots spatial representation of a microarray spot statistics. Currently
 the statistics can be either the log-ratios (default) or the intensities.
}

\usage{\method{plot}{MicroarrayData}(x, what, ..., style=NULL)}

\arguments{
 \item{what}{What to plot. Common formats are \code{"YvsX"} and
  \code{"spatialX"} where \code{X} and \code{Y} are variables in the
  object.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()}, \code{\link[aroma:plotSpatial.MicroarrayData]{*plotSpatial}()},
  \code{\link[graphics]{par}}, and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)
  rg <- as.RGData(ma)

  layout(matrix(1:6, nrow=2, byrow=TRUE))
  plot(raw)            # Plots R vs Rb. Default is "RvsRb" and slide 1.
  plot(raw, "GvsGb")   # Plots G vs Gb.
  plot(rg)             # Plots R vs G. Default is "RvsG" and slide 1.
  plot(rg, "GvsR")     # Plots G vs R.
  plot(ma, slide=4)    # Plots M vs A for slide 4. Default is "MvsA".
  plot(ma, "spatialM") # Plots spatial plot of M's.
}

\seealso{
  \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()} and \code{\link[aroma:plotSpatial.MicroarrayData]{*plotSpatial}()},
  which are called by \code{plot}.
  For highlightning and put labels on spots after plotting see
  \code{\link[aroma:highlight.MicroarrayData]{*highlight}()} and
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
