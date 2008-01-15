%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getForegroundSD.GenePixData}
\alias{getForegroundSD.GenePixData}
\alias{GenePixData.getForegroundSD}
\alias{getForegroundSD.GenePixData}
\alias{getForegroundSD,GenePixData-method}


\title{Gets the standard deviation of the foreground pixels}

\usage{\method{getForegroundSD}{GenePixData}(this, slides=NULL, channels=1:2, ...)}

\description{
 Gets the standard deviation of the foreground pixels.
}

\value{
  Returns a list of matrices that contain the estimated standard deviation
  of the pixels in the foreground region of the spots.
}

\details{
  The GenePix fields \code{"F532 SD"} (green) and \code{"F635 SD"} (red)
  are returned.
}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(gpr)
  sd <- getForegroundSD(gpr)
  raw$RSD <- sd$RSD; raw$GSD <- sd$GSD;
  sd <- getBackgroundSD(gpr)
  raw$RbSD <- sd$RbSD; raw$GbSD <- sd$GbSD;
  subplots(4)
  plot(raw, "RSDvsR", col="red")
  plot(raw, "RbSDvsRb", col="red")
  plot(raw, "GSDvsG", col="green")
  plot(raw, "GbSDvsGb", col="green")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
