%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotPrintorder.MicroarrayData}
\alias{plotPrintorder.MicroarrayData}
\alias{MicroarrayData.plotPrintorder}
\alias{plotPrintorder.MicroarrayData}
\alias{plotPrintorder,MicroarrayData-method}


\title{Plots the data as a time series in the order it was printed}

\description{
 Plots the data as a time series in the order it was printed.
}

\usage{\method{plotPrintorder}{MicroarrayData}(this, what, slides=NULL, include=NULL, exclude=NULL, cex="auto", col="auto", pch="auto", xlab="in print order", ylab=NULL, ylog=NULL, groupfcn=NULL, breakpoints=NULL, ..., f=200/nbrOfSpots(this), style=NULL)}

\arguments{
 \item{field}{The field to be plotted. Any field that can be retrieved
      by \code{extract}, is accepted.}
 \item{slides}{The slide(s) to be plotted. If \code{\link[base]{NULL}}, all slides are
      plotted.}
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
}

\details{
  If a vector of colors (\code{col}), symbol sizes (\code{cex}) or
  symbol types (\code{pch}) are specified, such vectors will be reordered
  according to the print order \emph{before} being applied.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  plotPrintorder(ma, "M")       # Printorder plot of log ratios before.
  normalizeWithinSlide(ma, "p") # Printtipwise lowess normalization.
  plotPrintorder(ma, "M")       # Printorder plot of log ratios after.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
 [1] Henrik Bengtsson, Plate Effects in cDNA microarray data,
     Matemathical Statistics, Centre for Matematical Sciences,
     Lund University, Sweden. Manuscript, 2002.
}

\seealso{
  \code{\link[aroma:plot.MicroarrayData]{*plot}()}.
  \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()}.
  \code{\link[aroma:highlight.MicroarrayData]{*highlight}()}.
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
