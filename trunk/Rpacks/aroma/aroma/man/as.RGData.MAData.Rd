%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.RGData.MAData}
\alias{as.RGData.MAData}
\alias{MAData.as.RGData}
\alias{as.RGData.MAData}
\alias{as.RGData,MAData-method}


\title{Transform MA format into RG format}

\description{
 Transform from log ratios and intensities to the red and green intensities.
 Background fields are set to zero.
}

\usage{\method{as.RGData}{MAData}(this, slides=NULL, ...)}

\arguments{
  \item{slides}{Subset of slides to be returned. If \code{\link[base]{NULL}}, all slides
  are returned.}
}

\value{
  Returns a \code{\link{RGData}} object.
}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  rg <- getSignal(raw, bgSubtract=TRUE)
  ma <- as.MAData(rg)
  rg2 <- as.RGData(ma)
  equal <- equals(rg, rg2)
}

\seealso{
  \code{\link{RGData.as.MAData}}.
  For more information see \code{\link{MAData}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
