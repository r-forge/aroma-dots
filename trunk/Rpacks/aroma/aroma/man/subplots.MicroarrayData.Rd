%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{subplots.MicroarrayData}
\alias{subplots.MicroarrayData}
\alias{MicroarrayData.subplots}
\alias{subplots.MicroarrayData}
\alias{subplots,MicroarrayData-method}


\title{Prepare for a grid of subplots based on the number of slides}

\description{
 Creates a grid of subplots based on the number of slides this dataset
 contains.
}

\usage{\method{subplots}{MicroarrayData}(this, ..., adjustMargins=TRUE, mar=c(2, 2, 1, 1))}

\arguments{
 \item{...}{Any arguments that \code{\link[R.graphics]{subplots.Device}} takes.}
}

\value{
  Returns a \code{\link[base]{matrix}} with the shape of the grids and which contains the
  subplot indices indicating the order that the subplots are plotted.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  subplots(ma)
  for (k in seq(ma))
    plot(ma, slide=k)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[graphics]{layout}}.
  \code{\link[R.graphics]{subplots.Device}}.
  To set the distances between the subplots see graphical parameter
  \code{mar} in \code{\link[graphics]{par}}.
}
\keyword{internal}
\keyword{methods}
