%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotMvsM.MAData}
\alias{plotMvsM.MAData}
\alias{MAData.plotMvsM}
\alias{plotMvsM.MAData}
\alias{plotMvsM,MAData-method}


\title{Plots the log-ratios for one slide versus another}

\description{
  Plots the log-ratios for one slide versus another.
}

\usage{\method{plotMvsM}{MAData}(this, slides=NULL, include=NULL, exclude=NULL, pch="auto", xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, ..., style=NULL)}

\arguments{
 \item{slides}{\code{\link[base]{vector}} of two slide indices to be plotted against
  each other.}
 \item{pch}{The plot symbol.}
 \item{xlim, ylim}{The visible range on the x and the y axis.}
 \item{xlab, ylab}{The labels of the x and the y axis.}
 \item{...}{Other arguments accepted by \code{\link[graphics]{plot}}.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
# Load some example data
SMA$loadData("mouse.data")
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
raw <- RawData(mouse.data, layout=layout)

# Get the signal
ma <- getSignal(raw, bgSubtract=TRUE)

subplots(4)

plotMvsM(ma, slides=c(1,2), xlim=c(-1,1)*3, pch=".")
points(ma, slides=c(3,4), col="red", pch=".")
points(ma, slides=c(5,6), col="blue", pch=".")
abline(a=0, b=1, lty=2)
title(main="Three different (M,M) pairs")


normalizeWithinSlide(ma, method="l")
normalizeAcrossSlides(ma)

plotMvsM(ma, slides=c(1,2), xlim=c(-1,1)*3, pch=".")
points(ma, slides=c(3,4), col="red", pch=".")
points(ma, slides=c(5,6), col="blue", pch=".")
abline(a=0, b=1, lty=2)
title(main="Same after lowess normalization")
}

\seealso{
  \code{\link[aroma:pointsMvsM.MAData]{*pointsMvsM}()}.
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
