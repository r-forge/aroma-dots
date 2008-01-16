%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% on Tue Jan 15 18:37:20 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getForegroundSE.SpotData}
\alias{getForegroundSE.SpotData}
\alias{SpotData.getForegroundSE}
\alias{getForegroundSE.SpotData}
\alias{getForegroundSE,SpotData-method}


\title{Gets the standard error of the foreground pixels}

\usage{\method{getForegroundSE}{SpotData}(this, ...)}

\description{
 Gets the standard error of the foreground pixels.
}

\value{
  Returns a list of matrices that contain the standard error
  of the pixels in the foreground region of the spots.
}

\details{
  The standard error returns the standard deviation divided by the area.
}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(gpr)
  sd <- getForegroundSE(gpr)
  raw$RSE <- sd$RSE; raw$GSE <- sd$GSE;
  subplots(4)
  plot(raw, "RSEvsR", col="red")
  plot(raw, "GSEvsG", col="green")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
