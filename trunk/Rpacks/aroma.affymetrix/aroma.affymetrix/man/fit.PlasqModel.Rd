%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PlasqModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fit.PlasqModel}
\alias{fit.PlasqModel}
\alias{PlasqModel.fit}
\alias{fit.PlasqModel}
\alias{fit,PlasqModel-method}


\title{Fits the model}

\description{
 Fits the model to all or to a subset of the units.
}

\usage{\method{fit}{PlasqModel}(this, units="remaining", ..., unitsPerChunk=moreUnits * 1e+05/length(getDataSet(this)), moreUnits=1, force=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns the indices of the units fitted, or \code{\link[base]{NULL}} if no units had to
 be fitted.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{PlasqModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
