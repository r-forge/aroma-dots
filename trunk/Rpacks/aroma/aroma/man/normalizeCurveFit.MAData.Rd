%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeCurveFit.MAData}
\alias{normalizeCurveFit.MAData}
\alias{MAData.normalizeCurveFit}
\alias{normalizeCurveFit.MAData}
\alias{normalizeCurveFit,MAData-method}


\title{Within-slide normalization that adjust log-ratios by estimating a smooth intensity-dependent curve in (A,M)}

\usage{\method{normalizeCurveFit}{MAData}(this, slides=NULL, groupBy=NULL, ..., .fits=NULL)}

\arguments{
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link{LayoutGroups}} specifying the
    groups of spots that are to normalized individually.
    If \code{\link[base]{NULL}}, all spots are normalized together.}
  \item{...}{Other arguments, such as \code{groupBy} and \code{method},
     passed to the curve estimator \code{\link[aroma:estimateMACurve.MAData]{*estimateMACurve}()}.}
}

\description{
  Within-slide normalization that adjust log-ratios by estimating a smooth intensity-dependent curve in (A,M).
}

\value{
  Returns itself invisibly.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
