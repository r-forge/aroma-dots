%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFlag.MicroarrayData}
\alias{getFlag.MicroarrayData}
\alias{MicroarrayData.getFlag}
\alias{getFlag.MicroarrayData}
\alias{getFlag,MicroarrayData-method}


\title{Gets the indices of the spots flagged with a specific flag}

\description{
  Gets the indices of the spots flagged with a specific flag.
}

\usage{\method{getFlag}{MicroarrayData}(this, flag, slide=NULL, ...)}

\arguments{
 \item{flag}{The name of the flag.}
 \item{slide}{Slide(s) to be considered. If \code{\link[base]{NULL}} all slides are
  considered.}
}

\value{
  Returns a \code{\link[base]{vector}} of \code{\link[base]{integer}}s.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)
  setFlag(ma, "|M|>1", abs(ma$M)>1)
  plot(ma)
  highlight(ma, "|M|>1", col="red")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:setFlag.MicroarrayData]{*setFlag}()},
  \code{\link[aroma:clearFlag.MicroarrayData]{*clearFlag}()},
  \code{\link[aroma:listFlags.MicroarrayData]{*listFlags}()},
  \code{\link[aroma:getInclude.MicroarrayData]{*getInclude}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
