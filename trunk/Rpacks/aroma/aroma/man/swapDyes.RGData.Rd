%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% on Tue Jan 15 18:36:16 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{swapDyes.RGData}
\alias{swapDyes.RGData}
\alias{RGData.swapDyes}
\alias{swapDyes.RGData}
\alias{swapDyes,RGData-method}


\title{Swap dyes of one or many slides}

\usage{\method{swapDyes}{RGData}(this, slides=NULL, ...)}

\description{
  Swap dyes of one or many slides.
}

\value{
  Return itself.
}

\arguments{
  \item{slides}{A \code{\link[base]{vector}} of slides to be dye swapped. If \code{\link[base]{NULL}}, all
  slides are considered.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)

  ma <- getSignal(raw)
  rg <- as.RGData(ma)

  # Dye swap every other slide.
  swapDyes(rg, slides=c(4,5,6))

  layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
  for (k in 1:6)
    plot(rg, slide=k)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{RGData}}.
}
\keyword{internal}
\keyword{methods}
