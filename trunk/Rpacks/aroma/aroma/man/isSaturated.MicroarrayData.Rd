%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:02 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{isSaturated.MicroarrayData}
\alias{isSaturated.MicroarrayData}
\alias{MicroarrayData.isSaturated}
\alias{isSaturated.MicroarrayData}
\alias{isSaturated,MicroarrayData-method}


\title{Check if some observations are saturated or not}

\description{
  Check if some observations are saturated or not.
}

\usage{\method{isSaturated}{MicroarrayData}(this, channels=NULL, index=NULL, slides=NULL, na.rm=TRUE, ...)}

\arguments{
  \item{channels}{\code{\link[base]{vector}} of channels to be checked. If \code{\link[base]{NULL}}, all are checked.}
  \item{index}{\code{\link[base]{vector}} of spot indicies to be checked. If \code{\link[base]{NULL}}, all are checked.}
  \item{slides}{\code{\link[base]{vector}} of slides to be checked. If \code{\link[base]{NULL}}, all are checked.}
  \item{na.rm}{If \code{\link[base:logical]{TRUE}}, observations for which it is unknown (==\code{\link[base]{NA}}) if
    it is saturated or not, are excluded.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:setSaturated.MicroarrayData]{*setSaturated}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
