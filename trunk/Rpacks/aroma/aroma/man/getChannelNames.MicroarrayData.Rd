%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getChannelNames.MicroarrayData}
\alias{getChannelNames.MicroarrayData}
\alias{MicroarrayData.getChannelNames}
\alias{getChannelNames.MicroarrayData}
\alias{getChannelNames,MicroarrayData-method}


\title{Gets the names of the channels}

\description{
 Gets the names of the channels.
}

\usage{\method{getChannelNames}{MicroarrayData}(this, ...)}

\arguments{
  \item{digitsOnly}{If \code{\link[base:logical]{TRUE}}, all non-digit characters are excluded from
    the channel names.}
  \item{settings}{Internal use only.}
}

\value{
 Returns an \code{\link[base]{character}} string \code{\link[base]{vector}}. An element with value \code{\link[base]{NA}} is a
 channel that had a zero-length name.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
