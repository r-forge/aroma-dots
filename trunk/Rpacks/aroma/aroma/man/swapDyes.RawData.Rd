%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% on Tue Jan 15 18:35:57 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{swapDyes.RawData}
\alias{swapDyes.RawData}
\alias{RawData.swapDyes}
\alias{swapDyes.RawData}
\alias{swapDyes,RawData-method}


\title{Dye swap one or many slides}

\usage{\method{swapDyes}{RawData}(this, slides=NULL, ...)}

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
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)

  # Dye swap every other slide.
  swapDyes(raw, slides=c(4,5,6))

  layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
  for (k in 1:6)
    plot(raw, "RvsG", slide=k)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{RawData}}.
}
\keyword{internal}
\keyword{methods}
