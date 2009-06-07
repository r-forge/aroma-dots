%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeLogRatioShift.MAData}
\alias{normalizeLogRatioShift.MAData}
\alias{MAData.normalizeLogRatioShift}
\alias{normalizeLogRatioShift.MAData}
\alias{normalizeLogRatioShift,MAData-method}


\title{Within-slide normalization that adjust the log-ratios shift}

\usage{\method{normalizeLogRatioShift}{MAData}(this, slides=NULL, groupBy=NULL, method=c("median", "mean"), ...)}

\arguments{
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link{LayoutGroups}} specifying the
    groups of spots that are to normalized individually.
    If \code{\link[base]{NULL}}, all spots are normalized together.}
  \item{method}{\code{\link[base]{character}} string specifying if the median or the mean
    should be used to estimate the shift of the log-ratios.}
}

\description{
  Within-slide normalization that adjust the log-ratios shift.
}

\value{
  Returns a \code{\link[base]{list}} structure containing information about the fit etc.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
