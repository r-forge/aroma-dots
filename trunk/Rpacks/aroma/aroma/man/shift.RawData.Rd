%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{shift.RawData}
\alias{shift.RawData}
\alias{RawData.shift}
\alias{shift.RawData}
\alias{shift,RawData-method}


\title{Shift the log-ratios, log-intensities or the raw signal}

\usage{\method{shift}{RawData}(this, M=NULL, A=NULL, R=NULL, G=NULL, Rb=NULL, Gb=NULL, slides=NULL, ...)}

\description{
  Shift the log-ratios, log-intensities or the raw signal.
}

\arguments{
  \item{M,A,R,G,Rb,Gb}{A \code{\link[base]{numeric}} or \code{\link[base]{function}} specifying the shift to be
   applied to the log-ratios, the log-intensities, the red, and/or the
   green foreground or background signals.
   If more than one of these are shifted at the same time, they are
   shifted in the order \code{M}, \code{A}, \code{R}, \code{G},
   \code{Rb} and \code{Gb}.
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
  For more information see \code{\link{RawData}}.
}
\keyword{internal}
\keyword{methods}
