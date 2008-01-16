%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{clearFlag.MicroarrayData}
\alias{clearFlag.MicroarrayData}
\alias{MicroarrayData.clearFlag}
\alias{clearFlag.MicroarrayData}
\alias{clearFlag,MicroarrayData-method}


\title{Clears the spots from the specified flag}

\description{
  Clears the spots from the specified flag. The argument \code{flag} can
  be a regular expression, which means that one can clear more than one
  flag at the time. All flags are cleared by \code{flag="*"}.
}

\usage{\method{clearFlag}{MicroarrayData}(this, flag, ...)}

\arguments{
 \item{flag}{The regular expression of which flags should be cleared.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getFlag.MicroarrayData]{*getFlag}()},
  \code{\link[aroma:setFlag.MicroarrayData]{*setFlag}()},
  \code{\link[aroma:listFlags.MicroarrayData]{*listFlags}()},
  \code{\link[aroma:getInclude.MicroarrayData]{*getInclude}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
