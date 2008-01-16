%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{text.MicroarrayData}
\alias{text.MicroarrayData}
\alias{MicroarrayData.text}
\alias{text.MicroarrayData}
\alias{text,MicroarrayData-method}


\title{Puts labels to the spots in last plot}

\description{
 Puts labels to the spots in last plot.
}

\usage{\method{text}{MicroarrayData}(x, include=NULL, exclude=NULL, labels=NULL, col="black", offset=0.2, pos=4, xpd=TRUE, ...)}

\arguments{
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded.
  If \code{\link[base]{NULL}} no spots are excluded.}
 \item{labels}{Labels of \emph{all} the data points. If \code{\link[base]{NULL}} first the labels of the belonging \link{Layout} object is retrieved. If there is no Layout object, the data points are labeled with indices of the plotted spots.}
 \item{col}{The color(s) to be used for the plotted spots, i.e. for the
  spots \emph{after} inclusion and exclusion.}
 \item{pos}{A position specifier for the text. Values of \code{1}, \code{2}, \code{3} and \code{4}, respectively indicate positions below, to the left of, above and to the right of the specified coordinates.}
 \item{offset}{When `pos' is specified, this value gives the offset of the label from the specified coordinate in fractions of a character width.}
 \item{xpd}{A logical value or \code{\link[base]{NA}}. If \code{\link[base:logical]{FALSE}}, all plotting is clipped to the plot region, if \code{\link[base:logical]{TRUE}}, all plotting is clipped to the figure region, and if \code{\link[base]{NA}}, all plotting is clipped to the device region.}
 \item{...}{Common arguments accepted by most text labeling functions.
  For more information see \code{\link[graphics]{text}} and \code{\link[graphics]{par}}.}
}

\examples{
  # Loads the file 'gpr123.gpr' located in the data directory:
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(gpr);
  ma <- getSignal(raw)

  # Highlight and label the first 50 genes
  idx <- 1:50
  plot(ma)

  # Highlight all spots in grid (1,2)
  highlight(ma, idx, col="purple")

  # Add bold faced labels using the labels specified by the layout
  text(ma, idx, font=2)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()}.
  \code{\link[aroma:plotSpatial.MicroarrayData]{*plotSpatial}()}.
  \code{\link[aroma:highlight.MicroarrayData]{*highlight}()}.
  \code{\link[graphics]{text}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
