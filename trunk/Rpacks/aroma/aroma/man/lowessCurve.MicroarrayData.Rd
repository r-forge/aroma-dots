%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{lowessCurve.MicroarrayData}
\alias{lowessCurve.MicroarrayData}
\alias{MicroarrayData.lowessCurve}
\alias{lowessCurve.MicroarrayData}
\alias{lowessCurve,MicroarrayData-method}


\title{Draws one or more lowess curves through the data in last plot}

\description{
 Draws one or more lowess curves through the data in last plot.
 Normally this function is only applicable to scatter plots.
}

\usage{\method{lowessCurve}{MicroarrayData}(this, include=NULL, exclude=NULL, col=NULL, f=0.3, ...)}

\arguments{
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded.
  If \code{\link[base]{NULL}} no spots are excluded.}
 \item{col}{The color(s) to be used for the plotted line(s). If \code{\link[base]{NULL}}
  the global default line color will be used.}
 \item{...}{Common arguments accepted by underlying plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  plot(ma)
  lowessCurve(ma, gridwise=TRUE)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:text.MicroarrayData]{*text}()}.
  \code{\link[graphics]{points}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
