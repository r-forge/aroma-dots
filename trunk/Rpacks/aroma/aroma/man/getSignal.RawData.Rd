%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSignal.RawData}
\alias{getSignal.RawData}
\alias{RawData.getSignal}
\alias{getSignal.RawData}
\alias{getSignal,RawData-method}


\title{Gets the red and the green signal}

\description{
  Gets the red and the green signal as a MAData object.
}

\usage{\method{getSignal}{RawData}(this, bgSubtract=TRUE, slides=NULL, bg.subtract=TRUE, ...)}

\arguments{
  \item{slides}{Subset of slides to be returned. If \code{\link[base]{NULL}}, all slides
  are returned.}
  \item{bgSubstract}{If \code{\link[base:logical]{TRUE}}, the background is subtracted from the
        foreground, before the transformation is performed.
        This argument was previously named \code{bg.subtract}, which still
        works, but is deprecated and will be removed in a future version.}
}

\value{
  Returns a \code{\link{MAData}} object.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=FALSE);
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getBackground.RawData]{*getBackground}()},
  \code{\link[aroma:getForeground.RawData]{*getForeground}()}.
  For more information see \code{\link{RawData}}.
}
\keyword{internal}
\keyword{methods}
