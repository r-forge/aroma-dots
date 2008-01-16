%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% on Tue Jan 15 18:36:16 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{shift.RGData}
\alias{shift.RGData}
\alias{RGData.shift}
\alias{shift.RGData}
\alias{shift,RGData-method}


\title{Shift the log-ratios, log-intensities or the raw signal}

\usage{\method{shift}{RGData}(this, M=NULL, A=NULL, R=NULL, G=NULL, slides=NULL, ...)}

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

\examples{\dontrun{See shift() for the MAData class for an example.}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{RGData}}.
}
\keyword{internal}
\keyword{methods}
