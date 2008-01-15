%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setLayout.MicroarrayData}
\alias{setLayout.MicroarrayData}
\alias{MicroarrayData.setLayout}
\alias{setLayout.MicroarrayData}
\alias{setLayout,MicroarrayData-method}


\title{Sets the layout}

\description{
  Sets the information about the layout structure of the microarray using
  a \link{Layout} object.
}

\usage{\method{setLayout}{MicroarrayData}(this, layout, ...)}

\arguments{
  \item{layout}{A new \link{Layout} object.}
}

\examples{
  SMA$loadData(c("mouse.data", "mouse.setup", "mouse.gnames"))
  raw <- RawData(mouse.data)
  setLayout(raw, as.Layout(mouse.setup, id=mouse.gnames))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
