%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% on Tue Jan 15 18:33:06 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeAffineShift.MAData}
\alias{normalizeAffineShift.MAData}
\alias{MAData.normalizeAffineShift}
\alias{normalizeAffineShift.MAData}
\alias{normalizeAffineShift,MAData-method}


\title{Within-slide normalization that adjust the channel bias}

\usage{\method{normalizeAffineShift}{MAData}(this, slides=NULL, groupBy=NULL, weights=NULL, method=c("Kerr", "Newton", "Bengtsson"), interval=c("safe", "auto"), ...)}

\arguments{
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{groupBy}{\code{\link[base]{character}} string or \code{\link{LayoutGroups}} specifying the
    groups of spots that are to normalized individually.
    If \code{\link[base]{NULL}}, all spots are normalized together.}
  \item{method}{\code{\link[base]{character}} string specifying which method to use.
    If \code{"Kerr"}, perpendicular translation normalization is applied.
    If \code{"Newton"}, parallel translation normalization is applied.
    If \code{"Bengtsson"}, semi-perpendicular translation normalization
    is applied.}
  \item{interval}{The interval of possible shifts.
    Either a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length 2 or a \code{\link[base]{character}} string.
    If \code{"safe"}, the interval will be set, depending on method,
    such that no signals are made non-positive.
    If \code{"auto"}, the interval is choosen automatically depending on
    method. For some methods, \code{"auto"} gives the same result
    as \code{"safe"}.}
  \item{...}{Currently not used.}
}

\description{
  Within-slide normalization that adjust the channel bias with the objective to minimize intensity dependency of
  the log-ratios.
}

\value{
  Returns a \code{\link[base]{list}} structure containing information about the fit etc.
}

\examples{
SMA$loadData("mouse.data")
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
raw <- RawData(mouse.data, layout=layout)
methods <- c("Kerr", "Bengtsson")
Alim <- c(6,16)
Mlim <- c(-1,1)*3
subplots(2*length(methods))
for (method in methods) {
  ma <- getSignal(raw, bgSubtract=TRUE)
  fit <- normalizeAffineShift(ma, method=method, slide=1)
  a <- fit$shift[1,1]
  msg <- sprintf("Optimal shift according to '", method, "' is a = \%1.f.\n", a)
  cat(msg)
  plot(ma, xlim=Alim, ylim=Mlim)
  lowessCurve(ma)
 
  x <- seq(-2*a,2*a, length=200)
  y <- sapply(x, FUN=fit$objectiveFunction, R=getR(ma),G=getG(ma))
  plot(x,y, pch=176, xlab="shift", ylab="objective");
  lines(x,y)
  abline(v=a)
} # for (method ...)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
