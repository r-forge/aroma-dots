%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Weights.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fromSaturatedSignals.SignalWeights}
\alias{fromSaturatedSignals.SignalWeights}
\alias{SignalWeights.fromSaturatedSignals}
\alias{fromSaturatedSignals.SignalWeights}
\alias{fromSaturatedSignals,SignalWeights-method}


\title{Weight function for scanner signals}

\synopsis{fromSaturatedSignals.SignalWeights(static, x, satRatio=0.75, maxSignal=2^16 - 1, ...)}
\usage{SignalWeights$fromSaturatedSignals(x, satRatio=0.75, maxSignal=2^16 - 1, ...)}

\description{
 Weight function for scanner signals.
}

\arguments{
 \item{x}{A \code{\link[base]{vector}} (or a \code{\link[base]{matrix}}) of scanner signals for which weights
   are calculated.}
 \item{satRatio}{A \code{\link[base]{double}} in [0,1] specifying the relative level (of the
   maximum signal) where saturations starts.}
 \item{maxSignal}{A \code{\link[base]{numeric}} specifying the maximum scanner signal.}
}

\value{
  Returns a \code{\link[base]{vector}} (or a \code{\link[base]{matrix}}) of weights in [0,1].
}

\details{
 All signals less than \code{satRatio*maxSignal} get weight one, all
 greater or equal to \code{maxSignal} get weight zero, and inbetween
 there is a linear decrease from one to zero.

 Default value of argument \code{satRatio} is chosen such that saturation
 "starts" about 50000 for signals on the range [0,65535], which was
 reported by [1].
}

\section{Missing values}{
  Missing-value signals (\code{\link[base]{NA}}), get zero weight.
}

\examples{
  x <- seq(from=0, to=10^5, by=100)
  plot(x, SignalWeights$fromSaturatedSignals(x), type="l")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
 [1] H. Lyng, A. Badiee, D.H. Svendsrud, E. Hovig, O. Myklebost, T. Stokke.
     Profound influence of microarray scanner characteristics on gene
     expression ratios: analysis and procedure for correction.
     BMC Genomics, 2004, 5:10.
}

\seealso{
  For more information see \code{\link{SignalWeights}}..
}
\keyword{internal}
\keyword{methods}
