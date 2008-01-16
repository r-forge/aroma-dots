%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{shift.MAData}
\alias{shift.MAData}
\alias{MAData.shift}
\alias{shift.MAData}
\alias{shift,MAData-method}


\title{Shift the log-ratios, log-intensities or the raw signal}

\usage{\method{shift}{MAData}(this, M=NULL, A=NULL, R=NULL, G=NULL, slides=NULL, ...)}

\description{
  Shift the log-ratios, log-intensities or the raw signal.
}

\arguments{
  \item{M,A,R,G}{A \code{\link[base]{numeric}} or \code{\link[base]{function}} specifying the shift to be
   applied to the log-ratios, the log-intensities, the red signals,
   and/or the green signals.
   If more than one of these are shifted at the same time, they are
   shifted in the order \code{M}, \code{A}, \code{R} and \code{G}.
   A \code{\link[base]{numeric}} specify the amount of shift.
   If a \code{\link[base]{function}}, e.g. \code{min()}, is used, then the amount of shift
   is the value returned by that function when all \emph{finite} values
   are passed to that function, e.g. \code{min(x[is.finite(x)])}. In
   other words, \code{\link[base]{NA}}'s etc are automatically taken care of.
  }
  \item{slides}{Slides to be shifted. If \code{\link[base]{NULL}}, all slides are shifted.}
}

\value{
  Returns nothing.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)

  subplots(4)
  xlim <- c(4,16); ylim <- c(-3,3);
  plot(ma, xlim=xlim, ylim=ylim)
  min1 <- function(x) { min(x)-1 }   # Shift to signal one (not zero!)
  shift(ma, R=min, G=min)
  plot(ma, xlim=xlim, ylim=ylim)
  shift(ma, M=median)
  plot(ma, xlim=xlim, ylim=ylim)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
