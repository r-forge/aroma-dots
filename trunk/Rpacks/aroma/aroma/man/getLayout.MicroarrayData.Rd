%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLayout.MicroarrayData}
\alias{getLayout.MicroarrayData}
\alias{MicroarrayData.getLayout}
\alias{getLayout.MicroarrayData}
\alias{getLayout,MicroarrayData-method}


\title{Gets the layout}

\usage{\method{getLayout}{MicroarrayData}(this, ...)}

\description{
  Gets the information about the layout structure of the microarray and
  return it as a \link{Layout} object.
}

\value{
  Returns a \code{\link{Layout}} object.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  print(getLayout(raw))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
