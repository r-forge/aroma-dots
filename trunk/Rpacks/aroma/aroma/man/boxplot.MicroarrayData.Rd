%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{boxplot.MicroarrayData}
\alias{boxplot.MicroarrayData}
\alias{MicroarrayData.boxplot}
\alias{boxplot.MicroarrayData}
\alias{boxplot,MicroarrayData-method}


\title{Plots a boxplot}

\description{
 Creates a box plot for any data over all grids or over all slides. If
 the argument \code{slide} is given, an boxplot over the print-tips will
 be generated, others an across-slides boxplot will be generated.
}

\usage{\method{boxplot}{MicroarrayData}(x, what, style=NULL, slides=NULL, include=NULL, exclude=NULL, groupBy="slide", gridwise=FALSE, names=NULL, xlab=NULL, ylab=NULL, log=NULL, las=3, col=NULL, ...)}

\arguments{
 \item{what}{What data variable to do a boxplot for.}
 \item{slides}{The slide(s) that should be used to generate this plot.
    If \code{\link[base]{NULL}} data from all slides are considered.}
  \item{include}{The indices of the spots that should be included.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are considered.
   If \code{\link[base]{NULL}} all spots are considered.}
  \item{exclude}{The indices of the spots that should be excluded.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are excluded.}
 \item{groupBy}{The data can be grouped by either \code{"slide"}
   or \code{"printtip"}.}
 \item{gridwise}{(Deprecated)
    If \code{\link[base:logical]{TRUE}} a within-slide (across-grid) boxplot
    is generated, otherwise an across-slide boxplot is generated.}
 \item{names}{A string vector of names for each grid/slide.
  If \code{\link[base]{NULL}} the grids will be named 1,2,3, etc.}
 \item{xlab, ylab}{Label on the x-axis (y-axis). If \code{\link[base]{NULL}}
  the label will be automatically set.}
 \item{las}{Rotation style of the axis labels. Default value is 3
  (always vertical). For more information see \code{\link[graphics]{par}}.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  boxplot(ma, groupBy="printtip") # Log-ratios for printtips on all slides
  boxplot(ma)                     # Log-ratios for all slides.
  boxplot(ma, slides=c(1,3,5))    # Log-ratios for slides 1, 3 and 5.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[graphics]{boxplot}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
