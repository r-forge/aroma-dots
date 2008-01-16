%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{hasLayout.MicroarrayData}
\alias{hasLayout.MicroarrayData}
\alias{MicroarrayData.hasLayout}
\alias{hasLayout.MicroarrayData}
\alias{hasLayout,MicroarrayData-method}


\title{Checks if the layout has been specified}

\usage{\method{hasLayout}{MicroarrayData}(this, ...)}

\description{
  Checks if the layout has been specified.
}

\value{
  Returns \code{\link[base:logical]{TRUE}} if a Layout object has been specified, otherwise \code{\link[base:logical]{FALSE}}.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  print(hasLayout(raw))  # TRUE
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
