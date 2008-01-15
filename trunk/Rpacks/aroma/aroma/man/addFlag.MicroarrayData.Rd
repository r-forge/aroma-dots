%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{addFlag.MicroarrayData}
\alias{addFlag.MicroarrayData}
\alias{MicroarrayData.addFlag}
\alias{addFlag.MicroarrayData}
\alias{addFlag,MicroarrayData-method}


\title{Flags or unflags the indices of some spots with a specific flag}

\description{
  Flags or unflags the indices of some spots with a specific flag.
  Note the difference from \code{\link[aroma:setFlag.MicroarrayData]{*setFlag}()} which sets
  flags regardless of previously flags set, whereas this add/removes to
  previously set flags.
}

\usage{\method{addFlag}{MicroarrayData}(this, flag, include=NULL, exclude=NULL, slide=NULL, include.op="and", exclude.op="or", ...)}

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

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getFlag.MicroarrayData]{*getFlag}()},
  \code{\link[aroma:setFlag.MicroarrayData]{*setFlag}()},
  \code{\link[aroma:clearFlag.MicroarrayData]{*clearFlag}()},
  \code{\link[aroma:listFlags.MicroarrayData]{*listFlags}()},
  \code{\link[aroma:getInclude.MicroarrayData]{*getInclude}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
