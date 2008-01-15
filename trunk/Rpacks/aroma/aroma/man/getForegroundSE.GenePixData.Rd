%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getForegroundSE.GenePixData}
\alias{getForegroundSE.GenePixData}
\alias{GenePixData.getForegroundSE}
\alias{getForegroundSE.GenePixData}
\alias{getForegroundSE,GenePixData-method}


\title{Gets the standard error of the foreground pixels}

\usage{\method{getForegroundSE}{GenePixData}(this, slides=NULL, channels=1:2, ...)}

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
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
