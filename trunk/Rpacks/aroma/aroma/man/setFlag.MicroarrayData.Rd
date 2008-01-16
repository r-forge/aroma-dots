%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setFlag.MicroarrayData}
\alias{setFlag.MicroarrayData}
\alias{MicroarrayData.setFlag}
\alias{setFlag.MicroarrayData}
\alias{setFlag,MicroarrayData-method}


\title{Sets the flags on some of the spots}

\description{
  Sets the flags on some of the spots.
  See also \code{\link[aroma:addFlag.MicroarrayData]{*addFlag}()} which adds/removes flags to
  already set flags.
}

\usage{\method{setFlag}{MicroarrayData}(this, flag, include=NULL, exclude=NULL, slide=NULL, include.op="and", exclude.op="or", ...)}

\arguments{
 \item{flag}{The name of the flag.}
 \item{include}{The indices of the spots that should be included.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are considered.
  If \code{\link[base]{NULL}} all spots are considered.}
 \item{exclude}{The indices of the spots that should be excluded.
  If it is instead a name of one or more flags, the spots which have been
  flagged with these flags are excluded.
  If \code{\link[base]{NULL}} no spots are excluded.}
 \item{slide}{Slide(s) to be considered. If \code{\link[base]{NULL}} all slides are
  considered.}
 \item{include.op}{If \code{"and"} a spot must be flagged by all the
  flags to be included, and if the value is \code{"or"} it enough if
  spot is flagged by only one of the flags.
  Note that this argument is only effective if the \code{include}
  argument(s) are flag names.}
 \item{exclude.op}{If \code{"and"} a spot must be flagged by all the
  flags to be excluded, and if the value is \code{"or"} it enough if
  spot is flagged by only one of the flags.
  Note that this argument is only effective if the \code{exclude}
  argument(s) are flag names.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)
  normalizeWithinSlide(ma, "s")
  tma <- as.TMAData(ma)
  setFlag(tma, "Top 5\% M's", top(abs(tma$M), 0.05))
  setFlag(tma, "Top 5\% T's", top(abs(tma$T), 0.05))
  setFlag(tma, "Top 5\% M's AND top 5\% T's", c("Top 5\% M's", "Top 5\% T's"))
  plot(tma)
  highlight(tma, "Top 5\% M's", col="red")
  highlight(tma, "Top 5\% T's", col="blue")
  highlight(tma, "Top 5\% M's AND top 5\% T's", col="purple")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:addFlag.MicroarrayData]{*addFlag}()},
  \code{\link[aroma:getFlag.MicroarrayData]{*getFlag}()},
  \code{\link[aroma:clearFlag.MicroarrayData]{*clearFlag}()},
  \code{\link[aroma:listFlags.MicroarrayData]{*listFlags}()},
  \code{\link[aroma:getInclude.MicroarrayData]{*getInclude}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
