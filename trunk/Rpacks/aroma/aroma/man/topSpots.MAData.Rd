%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{topSpots.MAData}
\alias{topSpots.MAData}
\alias{MAData.topSpots}
\alias{topSpots.MAData}
\alias{topSpots,MAData-method}


\title{Gets the top spots}

\usage{\method{topSpots}{MAData}(this, n=50, nbins=1, slide=1, include=NULL, exclude=NULL, limits=c("both", "upper", "lower"), X="A", Y="M", ...)}

\description{
 Gets (approximately) the \code{n} spots with the highest M values, if
 \code{n} greater or equal to one, and the 100*\code{n} percentage if
 \code{n} is less than one.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)
  normalizeWithinSlide(ma, "s")
  top100  <- topSpots(ma, 100);   # The top 100 spots
  top0.01 <- topSpots(ma, 0.01);  # The top 1 percent of the spots
}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
