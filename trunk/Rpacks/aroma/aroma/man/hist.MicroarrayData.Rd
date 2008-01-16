%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{hist.MicroarrayData}
\alias{hist.MicroarrayData}
\alias{MicroarrayData.hist}
\alias{hist.MicroarrayData}
\alias{hist,MicroarrayData-method}


\title{Plots a histogram}

\description{
 Creates a histogram plot for any data variable over all grids or over
 all slides. If the argument \code{slide} is given, an boxplot over the print-tips will
 be generated, others an across-slides boxplot will be generated.
}

\usage{\method{hist}{MicroarrayData}(x, what, slides=NULL, include=NULL, exclude=NULL, xlab=NULL, log=NULL, ...)}

\arguments{
  \item{what}{What data variable to do a histogram for.}
  \item{slides}{The slide(s) that should be used to generate this plot.
   If \code{\link[base]{NULL}} data from all slides are considered.}
  \item{include}{The indices of the spots that should be included.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are considered.
   If \code{\link[base]{NULL}} all spots are considered.}
  \item{exclude}{The indices of the spots that should be excluded.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are excluded. If \code{\link[base]{NULL}} no spots are excluded.}
  \item{xlab, ylab}{Label on the x-axis (y-axis). If \code{\link[base]{NULL}}
   the label will be automatically set.}
  \item{log}{The log base to be used for taking the logarithm of the data
    before generating the histogram. If \code{\link[base]{NULL}} the non-logarithm data
    is plotted.}
  \item{...}{Common arguments accepted by most plot functions.
    For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  hist(ma, what="M", slide=1)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[graphics]{boxplot}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
