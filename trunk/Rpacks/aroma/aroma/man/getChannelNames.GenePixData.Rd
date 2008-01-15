%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getChannelNames.GenePixData}
\alias{getChannelNames.GenePixData}
\alias{GenePixData.getChannelNames}
\alias{getChannelNames.GenePixData}
\alias{getChannelNames,GenePixData-method}


\title{Gets the names of the channels}

\description{
 Gets the names of the channels.
}

\usage{\method{getChannelNames}{GenePixData}(this, digitsOnly=TRUE, settings=NULL, ...)}

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
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
