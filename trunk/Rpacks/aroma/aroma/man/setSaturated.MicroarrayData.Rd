%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setSaturated.MicroarrayData}
\alias{setSaturated.MicroarrayData}
\alias{MicroarrayData.setSaturated}
\alias{setSaturated.MicroarrayData}
\alias{setSaturated,MicroarrayData-method}


\title{Set which observations in a given channel that are saturated or not}

\description{
  Set which observations in a given channel that are saturated or not (or unknown).
}

\usage{\method{setSaturated}{MicroarrayData}(this, channel, index=NULL, slides=NULL, status=TRUE, ...)}

\arguments{
  \item{channels}{\code{\link[base]{vector}} of channels to be set.
    If \code{\link[base]{NULL}}, all are considered.}
  \item{index}{\code{\link[base]{vector}} of spot indicies to be set.
    If \code{\link[base]{NULL}}, all are considered.}
  \item{slides}{\code{\link[base]{vector}} of slides to be set.
    If \code{\link[base]{NULL}}, all are considered.}
  \item{status}{A \code{\link[base]{matrix}} of size matching the \code{index} and
    \code{slides} arguments (if not the values will be looped over).
    If a value is \code{\link[base:logical]{TRUE}}, the corresponding observation is saturated.
    If \code{\link[base:logical]{FALSE}}, it is not and if \code{\link[base]{NA}}, it is unknown. }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:isSaturated.MicroarrayData]{*isSaturated}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
