%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% on Tue Jan 15 18:36:51 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.ScanAlyzeData}
\alias{getRawData.ScanAlyzeData}
\alias{ScanAlyzeData.getRawData}
\alias{getRawData.ScanAlyzeData}
\alias{getRawData,ScanAlyzeData-method}


\title{Gets the raw (foreground and background) intensites}

\description{
 Extracts the red and green spot intensitites (both foreground and background)
 and returns a \code{\link{RawData}} object.
}

\usage{\method{getRawData}{ScanAlyzeData}(this, slides=NULL, fg=c("auto"), bg=("auto"), ...)}

\arguments{
  \item{slides}{Specifying which slides to be extracted. If \code{\link[base]{NULL}},
    all slides are considered.}
  \item{fg}{If \code{"auto"}, the default foreground estimates according
    to \code{\link[aroma:getForeground.ScanAlyzeData]{*getForeground}()} is returned.}
  \item{bg}{If \code{"auto"}, the default foreground estimates according
    to \code{\link[aroma:getBackground.ScanAlyzeData]{*getBackground}()} is returned.}
}

\value{
  Returns a \code{\link{RawData}} object containing the specified slides.
}

\details{
  The R and Rb channels will come from the CH2* fields, and
  the G and Gb channels will come from the CH1* fields.
  To swap the channels just use dyeSwap().
}

\examples{
  sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))

  # Get the foreground and the background
  raw <- getRawData(sa)

  # The the background corrected data
  ma <- getSignal(raw, bgSubtract=FALSE)

  # Plot M vs A with a lowess line through the data points
  plot(ma, slide=1)
  lowessCurve(ma, lwd=2, gridwise=TRUE)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ScanAlyzeData}}.
}
\keyword{internal}
\keyword{methods}
