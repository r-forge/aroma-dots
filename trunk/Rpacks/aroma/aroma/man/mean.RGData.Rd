%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% on Tue Jan 15 18:36:16 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mean.RGData}
\alias{mean.RGData}
\alias{RGData.mean}
\alias{mean.RGData}
\alias{mean,RGData-method}


\title{Genewise Average Mean for channel R and G}

\description{
  Computes the genewise average mean across all slides in the
  \code{\link{RGData}} object.
}

\usage{\method{mean}{RGData}(x, inf.rm=TRUE, na.rm=TRUE, ...)}

\arguments{
 \item{inf.rm}{A \code{\link[base]{logical}} value indicating whether \code{\link[base:is.finite]{Inf}} values should be
        stripped before the computation proceeds.}
 \item{na.rm}{A \code{\link[base]{logical}} value indicating whether \code{\link[base]{NA}} values should be
        stripped before the computation proceeds.}
}

\value{Returns a new \code{\link{RGData}} object containing the result.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{RGData}}.
}
\keyword{internal}
\keyword{methods}
