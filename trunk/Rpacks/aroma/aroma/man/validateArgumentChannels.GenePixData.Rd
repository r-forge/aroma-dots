%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{validateArgumentChannels.GenePixData}
\alias{validateArgumentChannels.GenePixData}
\alias{GenePixData.validateArgumentChannels}
\alias{validateArgumentChannels.GenePixData}
\alias{validateArgumentChannels,GenePixData-method}


\title{Validates the argument 'channels'}

\description{
 Validates the argument 'channels'.
}

\usage{\method{validateArgumentChannels}{GenePixData}(this, channels, minLength=length(channels), maxLength=minLength, ...)}

\arguments{
  \item{channels}{Either an \code{\link[base]{integer}} or a \code{\link[base]{character}} string \code{\link[base]{vector}}.}
  \item{minLength, maxLength}{The minimum and the maximum length of the
    argument \code{channels}. These arguments are used internally.}
}

\value{
 Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} of the same length as \code{channels} of
 channel indicies that are known to be valid.
 If the argument was invalid an \code{\link[R.oo]{Exception}} is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
