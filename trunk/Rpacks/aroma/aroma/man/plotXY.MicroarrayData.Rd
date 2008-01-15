%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotXY.MicroarrayData}
\alias{plotXY.MicroarrayData}
\alias{MicroarrayData.plotXY}
\alias{plotXY.MicroarrayData}
\alias{plotXY,MicroarrayData-method}


\title{Plots a scatter plot of two fields}

\description{
 Plots a scatter plot of two fields.
 It is recommended to use the \code{plot} function instead of calling
 this method explicitly (see \code{\link[aroma:plot.MicroarrayData]{*plot}()}).
}

\usage{\method{plotXY}{MicroarrayData}(this, what, style=NULL, slide=NULL, include=NULL, exclude=NULL, xlog=NULL, ylog=NULL, xlab=NULL, ylab=NULL, cex=NULL, col="auto", pch="auto", gridwise=FALSE, f=0.3, lim=NULL, xlim=NULL, ylim=NULL, ...)}

\arguments{
 \item{what}{What to plot. Any two field that can be retrieved by
  \code{extract} are accepted.}
 \item{slide}{The slide to be plotted.}
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded. If \code{\link[base]{NULL}} no spots are excluded.}
 \item{gridwise}{If \code{\link[base:logical]{TRUE}} a lowess line for each grid will be drawn,
  otherwise only the global lowess line will be drawn.}
 \item{xlog, ylog}{The logarithmic base to be used to take the
  logarithm of the x and the y values, respectively. If \code{\link[base]{NULL}}, the
  non-logarithmic values are plotted.
  Note the difference from the definition in \code{\link[graphics]{par}}.}
 \item{xlab, ylab}{The labels on the x and the y axis, respectively. If
  \code{\link[base]{NULL}}, the fields default label will be used.}
 \item{cex}{The scale factor to be used. If \code{\link[base]{NULL}} the system default
  scale factor will be used.}
 \item{col}{The color(s) to be used for the plotted spots, i.e. for the
  spots \emph{after} inclusion and exclusion. If the value is
  \code{"redgreen"} a red to green palette is used.}
 \item{pch}{The dot style. Default value is \code{176}, which is a small
  circle.}
 \item{f}{The bandwidth for the lowess lines.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  plotXY(ma)                    # Plot M vs A (default)
  plotXY(ma, what=c("M","A"))   # Plot A vs M.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:plot.MicroarrayData]{*plot}()}.
  \code{\link[aroma:plotSpatial.MicroarrayData]{*plotSpatial}()}.
  \code{\link[aroma:highlight.MicroarrayData]{*highlight}()}.
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  \code{\link[R.graphics]{plotSymbols.Device}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
