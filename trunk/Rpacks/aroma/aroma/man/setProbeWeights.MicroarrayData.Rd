%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.WEIGHTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{setProbeWeights.MicroarrayData}
\alias{setProbeWeights.MicroarrayData}
\alias{MicroarrayData.setProbeWeights}
\alias{setProbeWeights.MicroarrayData}
\alias{setProbeWeights,MicroarrayData-method}


\title{Sets probe weights on one or several slides}

\usage{\method{setProbeWeights}{MicroarrayData}(this, weights, slides=NULL, .force=FALSE, ...)}

\description{
  Sets probe weights on one or several slides.
}

\arguments{
  \item{slides}{The slides which should be assign weights.
    If \code{\link[base]{NULL}}, all slides are included.}
  \item{weights}{A \code{\link[base]{matrix}} or a \code{\link[base]{vector}}. If a \code{\link[base]{matrix}} it must have
    number of rows equal to number of spots, and number of columns equal
    to the number of specified slides.
    If a \code{\link[base]{vector}}, it must be of length equal to the number of spots.
    If \code{\link[base]{NULL}}, the weights are removed.}
}

\value{
  Returns itself invisibly.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
