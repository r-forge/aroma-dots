%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.MAData.RGData}
\alias{as.MAData.RGData}
\alias{RGData.as.MAData}
\alias{as.MAData.RGData}
\alias{as.MAData,RGData-method}


\title{Transform from the red and green intensities into log ratios between them and the log product of them}

\description{
  Transform from the red and green intensities into log ratios between them and the log product of them.
}

\usage{\method{as.MAData}{RGData}(this, slides=NULL, ...)}

\arguments{
  \item{slides}{A \code{\link[base]{vector}} of \code{\link[base]{integer}}s indicating which slides to transformed. If \code{\link[base]{NULL}}, all slides are transformed.}
}

\value{
  Returns a \code{\link{MAData}} object.
}

\examples{
  SMA$loadData(c("mouse.data", "mouse.setup"))
  raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
  rg <- getSignal(raw)
  ma <- as.MAData(rg)
}

\seealso{
  \code{\link{MAData.as.RGData}}.
  For more information see \code{\link{RGData}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
