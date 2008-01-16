%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% on Tue Jan 15 18:33:06 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{estimateMACurve.MAData}
\alias{estimateMACurve.MAData}
\alias{MAData.estimateMACurve}
\alias{estimateMACurve.MAData}
\alias{estimateMACurve,MAData-method}


\title{Estimates a smooth intensity-dependent curve in (A,M)}

\usage{\method{estimateMACurve}{MAData}(this, slides=NULL, groupBy=NULL, weights=NULL, method=c("loess", "lowess", "spline", "robust.spline"), span=NULL, ...)}

\arguments{
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link{LayoutGroups}} specifying the
    groups of spots that are to normalized individually.
    If \code{\link[base]{NULL}}, all spots are normalized together.}
  \item{weights}{A \code{\link[base]{vector}} or \code{\link[base]{matrix}} of spot weights used when
    estimating the normalization function.}
  \item{method}{\code{\link[base]{character}} string specifying which method to use when
    fitting the intensity-dependent function.
    Supported methods:
     \code{"loess"} (better than lowess),
     \code{"lowess"} (classic; supports only zero-one weights),
     \code{"spline"} (more robust than lowess at lower and upper
                      intensities; supports only zero-one weights),
     \code{"robust.spline"} (better than spline).
  }
  \item{span}{A \code{\link[base]{double}} value specifying the bandwidth of the estimator used.}
  \item{...}{No used (allows other methods to pass additional garbage).}
}

\description{
  Estimates a smooth intensity-dependent curve in (A,M).
}

\value{
  Returns a \code{\link[base]{list}} where each element is a \code{\link[base]{list}} corresponding to a slide,
  which in turn contains estimates of each set of spots specified by
  the \code{groupBy} argument.
}

\section{Missing values}{
 The estimate will only be made based on complete finite observations.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Utilized by \code{\link[aroma:drawCurveFit.MAData]{*drawCurveFit}()} and \code{\link[aroma:normalizeCurveFit.MAData]{*normalizeCurveFit}()}.

  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
