%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{swapDyes.MAData}
\alias{swapDyes.MAData}
\alias{MAData.swapDyes}
\alias{swapDyes.MAData}
\alias{swapDyes,MAData-method}


\title{Dye swap one or many slides}

\usage{\method{swapDyes}{MAData}(this, slides=NULL, ...)}

\description{
  Dye swap one or many slides.
}

\value{
  Returns itself.
}

\arguments{
  \item{slides}{A \code{\link[base]{vector}} of slides to be dye swapped. If \code{\link[base]{NULL}}, all
  slides are considered.}
}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)

  ma <- getSignal(raw, bgSubtract=TRUE)
  rg <- as.RGData(ma)

  # Dye swap every other slide.
  swapDyes(ma, slides=c(4,5,6))

  layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
  for (k in 1:6)
    plot(ma, slide=k)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
